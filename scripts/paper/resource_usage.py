"""Low-overhead benchmark resource sampling; no recursive disk scans during runs."""
from pathlib import Path
import datetime, json, os, shutil, subprocess, threading, time

def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def descendants(root):
    found, pending = set(), [root]
    while pending:
        pid = pending.pop()
        if pid in found:
            continue
        found.add(pid)
        try:
            # Child processes may belong to any thread, not only the main thread.
            for task in (Path('/proc') / str(pid) / 'task').iterdir():
                pending.extend(int(x) for x in (task / 'children').read_text().split())
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            pass
    return found

def meminfo():
    return {key.rstrip(':'): int(value) * 1024 for key, value, *_ in
            (line.split() for line in Path('/proc/meminfo').read_text().splitlines())}

def scope_path(cmd):
    # systemd-run --scope detaches the workload from the CLI process ancestry.
    if not cmd or str(cmd[0]) != 'systemd-run':
        return None
    unit = next((str(x).split('=', 1)[1] for x in cmd if str(x).startswith('--unit=')), None)
    return Path('/sys/fs/cgroup/system.slice') / (unit + '.scope') if unit else None

def sample(pid, scope, scratch):
    pids = descendants(pid)
    cg = {}
    if scope is not None and scope.exists():
        for entry in scope.rglob('cgroup.procs'):
            try:
                pids.update(int(x) for x in entry.read_text().split())
            except FileNotFoundError:
                pass
        for key in ['memory.current', 'memory.peak', 'memory.swap.current']:
            try:
                cg[key] = int((scope / key).read_text())
            except (FileNotFoundError, ValueError):
                pass
    rss = swap = count = 0
    for child in pids:
        try:
            status = (Path('/proc') / str(child) / 'status').read_text()
            values = {line.split(':', 1)[0]: line.split(':', 1)[1].strip() for line in status.splitlines()}
            rss += int(values.get('VmRSS', '0 kB').split()[0]) * 1024
            swap += int(values.get('VmSwap', '0 kB').split()[0]) * 1024
            count += 1
        except (FileNotFoundError, ProcessLookupError):
            pass
    disk = shutil.disk_usage(scratch)
    memory = meminfo()
    stats = {}
    for line in Path('/proc/diskstats').read_text().splitlines():
        fields = line.split()
        if fields[2] == 'md0' or fields[2].startswith(('nvme', 'xvd')):
            stats[fields[2]] = [int(x) for x in fields[3:]]
    return dict(time=utc(), monotonic=time.monotonic(), process_count=count,
                process_tree_rss_bytes=rss, process_tree_swap_bytes=swap,
                cgroup=cg, host_mem_total_bytes=memory['MemTotal'],
                host_mem_available_bytes=memory['MemAvailable'],
                host_used_excluding_available_bytes=memory['MemTotal']-memory['MemAvailable'],
                host_cached_bytes=memory['Cached'], host_dirty_bytes=memory['Dirty'],
                disk_total_bytes=disk.total, disk_used_bytes=disk.used, disk_free_bytes=disk.free,
                diskstats=stats,
                cpu_stat=[line for line in Path('/proc/stat').read_text().splitlines() if line.startswith('cpu')],
                loadavg=Path('/proc/loadavg').read_text().strip(),
                cpu_pressure=Path('/proc/pressure/cpu').read_text().strip())

def run_measured(cmd, run_dir, *, scratch='/scratch', interval=5, **kwargs):
    run_dir = Path(run_dir)
    trace = run_dir / 'RESOURCE_USAGE.jsonl'
    summary_path = run_dir / 'RESOURCE_SUMMARY.json'
    assert not trace.exists() and not summary_path.exists(), 'Resource record already exists'
    baseline = sample(os.getpid(), None, scratch)
    proc = subprocess.Popen(cmd, **kwargs)
    stop = threading.Event()
    scope = scope_path(cmd)
    errors, peaks = [], {}
    first = None
    last = None
    n = 0
    with trace.open('x') as stream:
        def take_sample():
            nonlocal first, last, n
            try:
                row = sample(proc.pid, scope, scratch)
                stream.write(json.dumps(row, separators=(',', ':')) + '\n')
                stream.flush()
                first = first or row
                last = row
                n += 1
                for key in ['process_tree_rss_bytes', 'process_tree_swap_bytes',
                            'host_used_excluding_available_bytes', 'disk_used_bytes']:
                    peaks[key] = max(peaks.get(key, 0), row[key])
                for key, value in row['cgroup'].items():
                    peaks['cgroup_' + key] = max(peaks.get('cgroup_' + key, 0), value)
            except Exception as error:
                errors.append(repr(error))
        def sampler():
            take_sample()
            while not stop.wait(interval):
                take_sample()
        thread = threading.Thread(target=sampler, daemon=True)
        thread.start()
        try:
            rc = proc.wait()
        finally:
            stop.set()
            thread.join()
            take_sample()
    final_disk = shutil.disk_usage(scratch)
    data = dict(sampling_interval_seconds=interval, samples=n, errors=errors,
                started_utc=baseline['time'], finished_utc=utc(), sampled_peak=peaks,
                disk_used_before_bytes=baseline['disk_used_bytes'],
                disk_used_after_bytes=final_disk.used,
                disk_peak_growth_bytes=max(0, peaks.get('disk_used_bytes',0)-baseline['disk_used_bytes']),
                exit_code=rc,
                interpretation='RSS sums concurrent process resident sets and may double-count shared pages. '
                'GNU time max RSS remains in time.txt and is not a concurrent tree total. '
                'Disk values are filesystem-wide allocated space, including inputs and prior outputs; '
                'peak growth subtracts the pre-run baseline. Five-second sampled peaks can miss short spikes. '
                'Diskstats counters are cumulative per device/layer, not additive across RAID layers.')
    summary_path.write_text(json.dumps(data, indent=2) + '\n')
    return subprocess.CompletedProcess(cmd, rc)
