//---------------------------------------------------------------------------
// Copyright (C) 2021 Bo Xu <xubo123@big.ac.cn>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Command line options parser
//---------------------------------------------------------------------------

#include <cstring>
#include <iostream>
#include "optionparser.h"

namespace detail {

bool parseCommandLine(int argc, char* argv[], const Option optionList[], int optionListSize, Opt& opt, int& index, bool printErrorMsg, bool POSIXLY_CORRECT, bool reset) {
    static int optIndex = 1;
    static int charIndex = 1;
    static int nonOptIndex = 0;
    static bool stopOptionParsing = false;

    if (reset) {
        optIndex = 1;
        charIndex = 1;
        nonOptIndex = 0;
        stopOptionParsing = false;
    }

    if (optIndex >= argc)
        return false;

    bool ret = true;

    if (!stopOptionParsing && argv[optIndex][0] == '-' && argv[optIndex][1] != '-' && argv[optIndex][1] != 0) {
        int i;
        for (i = 0; i < optionListSize; i++) {
            if (optionList[i].name[0] == '-' && optionList[i].name[1] == argv[optIndex][charIndex] && optionList[i].name[2] == 0)
                break;
        }
        if (i < optionListSize) {
            if (optionList[i].hasArg) {
                if (argv[optIndex][charIndex + 1] == 0 && optIndex + 1 < argc) {
                    opt = { i, argv[++optIndex] };
                    index = ++optIndex;
                    charIndex = 1;
                    ret = true;
                }
                else if (argv[optIndex][charIndex + 1] != 0) {
                    opt = { i, argv[optIndex] + charIndex + 1 };
                    index = ++optIndex;
                    charIndex = 1;
                    ret = true;
                }
                else {
                    if (printErrorMsg)
                        std::cerr << "Error: option '" << optionList[i].name << "' requires an argument.\n";
                    ret = false;
                }
            }
            else {
                opt = { i, nullptr };
                if (argv[optIndex][++charIndex] == 0) {
                    optIndex++;
                    charIndex = 1;
                }
                index = optIndex;
                ret = true;
            }
        }
        else {
            if (printErrorMsg)
                std::cerr << "Error: unrecognized option '" << '-' << argv[optIndex][charIndex] << "'\n";
            ret = false;
        }
    }
    else if (!stopOptionParsing && argv[optIndex][0] == '-' && argv[optIndex][1] == '-' && argv[optIndex][2] != '=' && argv[optIndex][2] != 0) {
        std::size_t optSize = std::strlen(argv[optIndex]);
        std::size_t iOptSize;
        int i;
        for (i = 0; i < optionListSize; i++) {
            iOptSize = std::strlen(optionList[i].name);
            if (optSize >= iOptSize && std::memcmp(argv[optIndex], optionList[i].name, iOptSize) == 0 && (argv[optIndex][iOptSize] == 0 || argv[optIndex][iOptSize] == '='))
                break;
        }
        if (i < optionListSize) {
            if (optionList[i].hasArg) {
                if (argv[optIndex][iOptSize] == 0 && optIndex + 1 < argc) {
                    opt = { i, argv[++optIndex] };
                    index = ++optIndex;
                    ret = true;
                }
                else if (argv[optIndex][iOptSize] == '=' && argv[optIndex][iOptSize + 1] != 0) {
                    opt = { i, argv[optIndex] + iOptSize + 1 };
                    index = ++optIndex;
                    ret = true;
                }
                else {
                    if (printErrorMsg)
                        std::cerr << "Error: option '" << optionList[i].name << "' requires an argument.\n";
                    ret = false;
                }
            }
            else {
                if (argv[optIndex][iOptSize] == 0) {
                    opt = { i, nullptr };
                    index = ++optIndex;
                    ret = true;
                }
                else {
                    if (printErrorMsg)
                        std::cerr << "Error: option '" << optionList[i].name << "' doesn't allow an argument.\n";
                    ret = false;
                }
            }
        }
        else {
            if (printErrorMsg)
                std::cerr << "Error: unrecognized option '" << argv[optIndex] << "'\n";
            ret = false;
        }
    }
    else {
        if (!stopOptionParsing) {
            if (argv[optIndex][0] == '-' && argv[optIndex][1] == '-' && argv[optIndex][2] == 0) {
                stopOptionParsing = true;
                index = ++optIndex;
                if (optIndex >= argc)
                    ret = false;
            }
            else if (POSIXLY_CORRECT)
                stopOptionParsing = true;
        }

        if (ret) {
            int i;
            for (i = nonOptIndex; i < optionListSize; i++) {
                if (optionList[i].name[0] == 0)
                    break;
            }
            if (i < optionListSize) {
                opt = { i, argv[optIndex] };
                index = ++optIndex;
                nonOptIndex = i + 1;
                ret = true;
            }
            else {
                if (printErrorMsg)
                    std::cerr << "Error: unrecognized argument '" << argv[optIndex] << "'\n";
                ret = false;
            }
        }
    }

    return ret;
}

};

int parseCommandLine(int argc, char* argv[], const Option optionList[], int optionListSize, Opt* opt, int* index, int printErrorMsg, int POSIXLY_CORRECT, int reset) {
    return detail::parseCommandLine(argc, argv, optionList, optionListSize, *opt, *index, !!printErrorMsg, !!POSIXLY_CORRECT, !!reset);
}
