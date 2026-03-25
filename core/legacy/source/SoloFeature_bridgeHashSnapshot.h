#ifndef H_SoloFeature_bridgeHashSnapshot
#define H_SoloFeature_bridgeHashSnapshot

class Parameters;

namespace solo_bridge_hash_snapshot {

// True when replay should skip the mapChunk read loop (requires SNAP_IN + skip-reads env).
bool replaySkipReadsEnabled(const Parameters &P);
bool stopAfterCountEnabled(const Parameters &P);

} // namespace solo_bridge_hash_snapshot

#endif
