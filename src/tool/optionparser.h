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

#ifndef OPTIONPARSER_H
#define OPTIONPARSER_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    const char* name;
    int hasArg;
    int val;
} Option;

typedef struct {
    int index;
    const char* val;
} Opt;

int parseCommandLine(int argc, char* argv[], const Option optionList[], int optionListSize, Opt* opt, int* index, int printErrorMsg, int POSIXLY_CORRECT, int reset);

#ifdef __cplusplus
}
#endif

#endif
