#include "ReadAlign.h"

#include "ErrorWarning.h"
#include "input/CbqStarAdapter.h"

#include <iostream>
#include <sstream>

int ReadAlign::loadCbqReadView(const star::input::CbqReadView& view) {
    star::input::CbqStarAdapterOptions options;
    options.read_nends = P.readNends;
    options.out_sam_read_id_number = P.outSAMreadIDnumber;
    options.out_qs_conversion_add = P.outQSconversionAdd;
    options.trim_cutadapt_enabled = (P.trimCutadapt == "Yes");
    options.preserve_read_name_extra = false;

    star::input::CbqStarReadBuffers buffers;
    buffers.read_name_mates = readNameMates;
    buffers.read0 = Read0;
    buffers.read1 = Read1;
    buffers.qual0 = Qual0;
    buffers.read_name_extra = &readNameExtra;
    buffers.read_length = readLength;
    buffers.read_length_original = readLengthOriginal;
    buffers.i_read_all = &iReadAll;
    buffers.read_files_index = &readFilesIndex;
    buffers.read_filter = &readFilter;
    buffers.read_file_type = &readFileType;

    std::string error;
    if (!star::input::load_cbq_read_view_into_star_mates(view, options, &buffers, &clipMates, &error)) {
        std::ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR in CBQ input adapter\n";
        errOut << error << "\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    }

    return readFileType;
}
