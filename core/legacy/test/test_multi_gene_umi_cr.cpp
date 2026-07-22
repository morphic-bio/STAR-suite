#include "MultiGeneUmiCr.h"

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace mg = multi_gene_umi_cr;

int main()
{
    mg::Result result = mg::resolve({{7, 5, 2}});
    assert(result.accepted && result.gene == 7);

    result = mg::resolve({{3, 10, 4}, {8, 9, 4}});
    assert(result.accepted && result.gene == 3);

    result = mg::resolve({{3, 10, 4}, {8, 10, 3}});
    assert(!result.accepted && result.reason == "corrected_count_tie");

    // Gene 3 wins after correction, but gene 8 has more reads whose original
    // UMI was already the corrected UMI. STAR rejects this molecule.
    result = mg::resolve({{3, 10, 2}, {8, 9, 3}});
    assert(!result.accepted && result.reason == "original_umi_dominance");

    result = mg::resolve({{3, 10, 3}, {8, 9, 3}});
    assert(result.accepted && result.gene == 3);

    bool threw = false;
    try { mg::resolve({{1, 1, 1}, {1, 2, 1}}); }
    catch (const std::invalid_argument &) { threw = true; }
    assert(threw);
    std::cout << "MultiGeneUMI_CR helper tests passed\n";
    return 0;
}
