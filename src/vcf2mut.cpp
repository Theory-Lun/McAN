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
// Read mutations from VCF file and/or convert VCF into mutation format file
//     The mutation format:
//         <SampleName>\t<AccessionID>\t[<POS>(<VariantType>:<REF>-><ALT>)[;<POS>(<VariantType>:<REF>-><ALT>)]...]
//---------------------------------------------------------------------------

#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef NDEBUG
#define ASSERT(e) ((void)0)
#else
#define ASSERT(e) ((!!(e)) ? (void)0 : [] { assert(!#e); }())
#endif

typedef int(*Callback)(const char*, size_t, size_t, size_t, void*);

extern "C" {
    int vcf2mut(const char* vcfName, const char* mutName, Callback lineHandler, void* arg, const char* missingValue);
}

namespace detail {

template<typename T>
inline constexpr const T& min(const T& x, const T& y) noexcept {
    return ((y < x) ? y : x);
}

template<typename T>
inline constexpr const T& max(const T& x, const T& y) noexcept {
    return ((x < y) ? y : x);
}

template<typename T>
inline constexpr T abs(const T& x) noexcept {
    return ((x < 0) ? -x : x);
}

inline constexpr bool charCaseInsCmp(char x, char y) noexcept {
    return (x == y || abs<char>(x - y) == 'a' - 'A');
}

inline constexpr bool isalpha(char x) noexcept {
    return ((x >= 'A' && x <= 'Z') || (x >= 'a' && x <= 'z'));
}

inline constexpr char code(char x, char y) noexcept {
    return (ASSERT(isalpha(x) && isalpha(y)),
        "KWRYSM"[(x - ((x < 'a') ? 'A' : 'a') + y - ((y < 'a') ? 'A' : 'a') + 3) % 7]);
}

class Substring {
public:
    typedef char value_type;
    typedef std::char_traits<value_type> traits_type;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef const value_type* const_pointer;
    typedef const value_type& const_reference;

    static constexpr size_type npos = size_type(-1);

    constexpr Substring() noexcept : _size(0), _pstr(nullptr) {}
    constexpr Substring(const char* str, size_type count) noexcept : _size(count), _pstr(str) {}
    template<size_type N>
    constexpr Substring(const char(&str)[N]) noexcept : _size(N - 1), _pstr(str) {}
    constexpr Substring(std::nullptr_t) = delete;
    Substring(const std::string& str) : _size(str.size()), _pstr(str.data()) {}

    operator std::string() const {
        return std::string(this->_pstr, this->_size);
    }

    constexpr const_reference operator [] (size_type pos) const noexcept {
        return (ASSERT(this->_pstr && this->_size > pos), this->_pstr[pos]);
    }
    constexpr const_reference front() const noexcept {
        return (ASSERT(this->_pstr && this->_size > 0), this->_pstr[0]);
    }
    constexpr const_reference back() const noexcept {
        return (ASSERT(this->_pstr && this->_size > 0), this->_pstr[this->_size - 1]);
    }
    constexpr const_pointer data() const noexcept { return this->_pstr; }
    constexpr size_type size() const noexcept { return this->_size; }
    constexpr size_type length() const noexcept { return this->_size; }
    constexpr bool empty() const noexcept { return (this->_size == 0); }

    constexpr Substring substr(size_type pos = 0, size_type count = npos) const noexcept {
        return ((this->_pstr) ?
            Substring(this->_pstr + min(this->_size, pos), min(count, max(this->_size, pos) - pos)) :
            Substring());
    }

    size_type find_first_of(const char* str, size_type pos, size_type count) const noexcept {
        if (str && count && this->_pstr) {
            for (; pos < this->_size; pos++)
                if (traits_type::find(str, count, this->_pstr[pos]))
                    return pos;
        }
        return npos;
    }
    size_type find_first_of(char ch, size_type pos = 0) const noexcept {
        const char* p = (this->_pstr && this->_size) ?
            traits_type::find(this->_pstr + min(this->_size, pos), max(this->_size, pos) - pos, ch) : nullptr;
        return ((p) ? p - this->_pstr : npos);
    }
    size_type find_first_of(Substring str, size_type pos = 0) const noexcept {
        return this->find_first_of(str.data(), pos, str.size());
    }

    template<typename Container>
    void split(Container& cont, size_type n, char delimitor = ' ', bool consecutiveDelimitersAsOne = false) const;

    template<typename Container>
    void split(Container& cont, char delimitor = ' ', bool consecutiveDelimitersAsOne = false) const {
        split(cont, npos, delimitor, consecutiveDelimitersAsOne);
    }

    void print(std::ostream& out) const {
        out.write(this->_pstr, this->_size);
    }

private:
    size_type _size;
    const_pointer _pstr;
};

template<typename Container>
void Substring::split(Container& cont, size_type n, char delimitor, bool consecutiveDelimitersAsOne) const {
    if (n == 0)
        return;

    if (!this->_pstr) {
        cont.emplace_back(Substring());
        return;
    }

    size_type s = 0;
    size_type i = 0;

    for (; i < this->_size; i++) {
        if (this->_pstr[i] == delimitor) {
            cont.emplace_back(this->_pstr + s, i - s);
            if (--n == 0)
                return;
            if (consecutiveDelimitersAsOne)
                for (; i + 1 < this->_size && this->_pstr[i + 1] == delimitor; i++);
            s = i + 1;
        }
    }
    cont.emplace_back(this->_pstr + s, i - s);
}

inline bool operator == (const Substring& left, const Substring& right) noexcept {
    return (left.size() == right.size() &&
        Substring::traits_type::compare(left.data(), right.data(), left.size()) == 0);
}

inline bool operator != (const Substring& left, const Substring& right) noexcept {
    return !(left == right);
}

inline std::ostream& operator << (std::ostream& out, const Substring& str) {
    str.print(out);
    return out;
}

class Variant {
public:
    enum Type { None = 0, Deletion, Insertion, SNP, Indel };
    static const char* str[];

    Variant() = default;
    Variant(unsigned long pos, Substring ref, Substring alt);

    Type type() const noexcept { return this->_type; }
    void print(std::ostream& out) const {
        out << this->_pos << '(' << str[this->_type] << ':' << this->_ref << "->" << this->_alt << ')';
    }

private:
    unsigned long _pos;
    Type _type;
    std::string _ref;
    std::string _alt;
};

const char* Variant::str[] = { "None", "Deletion", "Insertion", "SNP", "Indel" };

Variant::Variant(unsigned long pos, Substring ref, Substring alt) {
    std::size_t i = 0;
    auto refLen = ref.size();
    auto altLen = alt.size();
    auto minLen = min(refLen, altLen);

    for (; i < minLen; i++)
        if (!charCaseInsCmp(ref[i], alt[i]))
            break;

    this->_pos = pos + i;
    auto newRef = ref.substr(i);
    auto newAlt = alt.substr(i);

    if (i == minLen) {
        if (minLen < refLen) {
            this->_type = Deletion;
            newAlt = Substring("-");
        }
        else if (minLen < altLen) {
            this->_type = Insertion;
            this->_pos--;
            newRef = Substring("-");
        }
        else
            this->_type = None;
    }
    else if (i == minLen - 1 && refLen == altLen)
        this->_type = SNP;
    else if (i == 0 && minLen == 1 && pos == 1 && charCaseInsCmp(ref.back(), alt.back())) {
        if (minLen < refLen) {
            this->_type = Deletion;
            newRef = ref.substr(0, refLen - 1);
            newAlt = Substring("-");
        }
        else if (minLen < altLen) {
            this->_type = Insertion;
            this->_pos--;
            newRef = Substring("-");
            newAlt = alt.substr(0, altLen - 1);
        }
    }
    else
        this->_type = Indel;

    this->_ref = newRef;
    this->_alt = newAlt;
}

inline std::ostream& operator << (std::ostream& out, const Variant& variant) {
    variant.print(out);
    return out;
}

typedef std::vector<Variant> Variants;

inline void printVariants(std::ostream& out, const Variants& variants) {
    bool first = true;
    for (const auto& variant : variants) {
        if (variant.type() == Variant::None)
            continue;
        if (!first)
            out << ';';
        else
            first = false;
        out << variant;
    }
}

inline std::ostream& operator << (std::ostream& out, const Variants& variants) {
    printVariants(out, variants);
    return out;
}

enum VCFField { CHROM = 0, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE };

struct VCFContent {
    std::string header;
    std::vector<Substring> headerFields;
    std::vector<Variants> sampleVariants;
};

bool parseVCF(std::istream& in, VCFContent& vcf) {
    std::string line;
    std::vector<Substring> fields;
    std::vector<Substring> alleles;
    std::size_t sampleSize;
    std::size_t recordSize = 0;
    bool headerHasRead = false;

    while (std::getline(in, line)) {
        if (line.size() >= 1 && line[0] == '#') {
            if (headerHasRead) {
                std::cerr << "\nError: invalid VCF file.\n";
                return false;
            }

            if (line.size() >= 2 && line[1] == '#')
                continue;

            vcf.header = line;
            Substring(vcf.header).split(vcf.headerFields, '\t');
            if (vcf.headerFields.size() <= SAMPLE) {
                std::cerr << "\nWarning: nothing to do.\n";
                return true;
            }

            sampleSize = vcf.headerFields.size() - SAMPLE;
            vcf.sampleVariants.resize(sampleSize);
            headerHasRead = true;
            continue;
        }

        fields.clear();
        Substring(line).split(fields, '\t');
        if (fields.size() != vcf.headerFields.size()) {
            std::cerr << "\nError: invalid VCF file.\n";
            return false;
        }

        if (fields[ALT] == "*" || fields[ALT] == ".")
            continue;

        alleles.clear();
        alleles.push_back(fields[REF]);
        fields[ALT].split(alleles, ',');

        if (fields[FORMAT].substr(0, 2) != "GT" || (fields[FORMAT].size() > 2 && fields[FORMAT][2] != ':'))
            continue;

        char* pEnd = nullptr;
        auto variantPos = std::strtoul(fields[POS].data(), &pEnd, 10);
        if (fields[POS].data() + fields[POS].size() != pEnd) {
            std::cerr << "Warning: invalid POS field in record " << recordSize + 1 << ". Skipped.\n";
            continue;
        }

        for (std::size_t i = SAMPLE; i < fields.size(); i++) {
            auto pos1 = fields[i].find_first_of(':');
            auto pos2 = fields[i].substr(0, pos1).find_first_of("/|");

            if (pos2 == Substring::npos) {
                auto allele = fields[i].substr(0, pos1);

                if (allele == "." || allele == "0")
                    continue;

                pEnd = nullptr;
                auto idx = static_cast<std::size_t>(std::strtoul(allele.data(), &pEnd, 10));
                if (allele.data() + allele.size() != pEnd || idx >= alleles.size()) {
                    std::cerr << "Warning: invalid genotype in record " << recordSize + 1 << ". Skipped.\n";
                    continue;
                }

                vcf.sampleVariants[i - SAMPLE].emplace_back(variantPos, alleles[0], alleles[idx]);
            }
            else {
                auto pos3 = fields[i].substr(0, pos1).find_first_of("/|", pos2 + 1);
                if (pos3 != Substring::npos) {
                    std::cerr << "Warning: polyploidy genotype in record " << recordSize + 1 << ". Haplotypes other than the first two are ignored.\n";
                    pos1 = pos3;
                }

                auto allele1 = fields[i].substr(0, pos2);
                auto allele2 = fields[i].substr(pos2 + 1, pos1 - (pos2 + 1));

                if ((allele1 == "." && allele2 == ".") || (allele1 == "0" && allele2 == "0"))
                    continue;

                std::size_t idx1 = 0;
                if (allele1 != "." && allele1 != "0") {
                    pEnd = nullptr;
                    idx1 = static_cast<std::size_t>(std::strtoul(allele1.data(), &pEnd, 10));
                    if (allele1.data() + allele1.size() != pEnd || idx1 >= alleles.size()) {
                        std::cerr << "Warning: invalid genotype in record " << recordSize + 1 << ". Skipped.\n";
                        continue;
                    }
                }

                std::size_t idx2 = 0;
                if (allele2 != "." && allele2 != "0") {
                    pEnd = nullptr;
                    idx2 = static_cast<std::size_t>(std::strtoul(allele2.data(), &pEnd, 10));
                    if (allele2.data() + allele2.size() != pEnd || idx2 >= alleles.size()) {
                        std::cerr << "Warning: invalid genotype in record " << recordSize + 1 << ". Skipped.\n";
                        continue;
                    }
                }

                if (idx1 == idx2)
                    vcf.sampleVariants[i - SAMPLE].emplace_back(variantPos, alleles[0], alleles[idx1]);
                else if (alleles[0].size() == 1 && alleles[idx1].size() == 1 && alleles[idx2].size() == 1)
                    vcf.sampleVariants[i - SAMPLE].emplace_back(variantPos, alleles[0], std::string(1, code(alleles[idx1][0], alleles[idx2][0])));
                else {
                    if (idx1 != 0 && idx2 != 0)
                        std::cerr << "Warning: genotype consisting of two distinct ATL alleles in record " << recordSize + 1 << ". The first haplotype is ignored.\n";
                    vcf.sampleVariants[i - SAMPLE].emplace_back(variantPos, alleles[0], alleles[(idx2 == 0) ? idx1 : idx2]);
                }
            }
        }

        recordSize++;
    }

    return true;
}

}

int vcf2mut(const char* vcfName, const char* mutName, Callback lineHandler, void* arg, const char* missingValue) {
    std::ifstream input;
    if (vcfName[0] != '-' || vcfName[1] != 0) {
        input.open(vcfName);
        if (!input.good()) {
            std::cerr << "\nError: failed to open file: " << vcfName << '\n';
            return -1;
        }
    }
    std::istream& in = (vcfName[0] != '-' || vcfName[1] != 0) ? input : std::cin;

    detail::VCFContent vcf;
    int ret = 0;

    if (!detail::parseVCF(in, vcf))
        return -1;

    auto sampleSize = vcf.sampleVariants.size();

    if (mutName) {
        std::ofstream out;
        out.open(mutName);
        if (!out.good()) {
            std::cerr << "\nError: failed to create file: " << mutName << '\n';
            return -1;
        }

        for (std::size_t i = 0; i < sampleSize; i++) {
            if (missingValue)
                out << missingValue;
            out << '\t' << vcf.headerFields[i + detail::SAMPLE] << '\t' << vcf.sampleVariants[i] << '\n';
        }
    }

    if (lineHandler) {
        std::ostringstream out;
        for (std::size_t i = 0; i < sampleSize; i++) {
            out.str("");
            out.clear();
            if (missingValue)
                out << missingValue;
            out << '\t' << vcf.headerFields[i + detail::SAMPLE] << '\t' << vcf.sampleVariants[i];

            vcf.sampleVariants[i].clear();
            vcf.sampleVariants[i].shrink_to_fit();

            auto line = out.str();
            ret = lineHandler(line.c_str(), line.size(), i, sampleSize, arg);
            if (ret)
                break;
        }
    }

    return ret;
}
