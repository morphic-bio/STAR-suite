#include "FlexDebugCounters.h"

#ifdef DEBUG_FLEX_COUNTERS

// Thread-local counter instance
thread_local FlexDebugCounters g_flexCounters;

// Global registry
std::vector<FlexDebugCounters*> g_allFlexCounters;
std::mutex g_flexCountersMutex;

void flexCountersRegisterThread() {
    std::lock_guard<std::mutex> lock(g_flexCountersMutex);
    // Check if this thread's counters are already registered
    for (auto* p : g_allFlexCounters) {
        if (p == &g_flexCounters) return;
    }
    g_allFlexCounters.push_back(&g_flexCounters);
}

FlexDebugCounters flexCountersSumAll() {
    FlexDebugCounters total;
    std::lock_guard<std::mutex> lock(g_flexCountersMutex);
    for (auto* c : g_allFlexCounters) {
        c->addTo(total);
    }
    return total;
}

#endif // DEBUG_FLEX_COUNTERS
