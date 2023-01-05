# McAN


## Introduction

McAN is a tool for haplotype network construction for sequences. The web service of McAN will be found here https://ngdc.cncb.ac.cn/ncov/online/tool/haplotype. The haplotype network of SARS-CoV-2 constructed by McAN will be found here https://ngdc.cncb.ac.cn/ncov/haplotype/.

## Compiling

McAN can be build on most Linux systems. It depends on two third-party libraries:
- cJSON: https://github.com/DaveGamble/cJSON.git
- hashmap: https://github.com/DavidLeeds/hashmap.git

Please install
- cmake ( >= 3.7 )
for compiling and running McAN. And then compile McAN by doing follows:
```shell
git clone https://github.com/Theory-Lun/McAN.git
cd McAN
mkdir build
cd build
cmake ..
make
```

## Usage and Examples

For testing McAN, please doing this:
```shell
cd .../McAN/example
mkdir out
```
First, filter sites and generate a sitemask file:
```shell
../build/McAN siteMask \
  --mutation data/mut \
  --meta data/meta \
  --out out/siteMask \
  --minfreq 0.01
```
Then, construct haplotype network:
```shell
../build/McAN \
  --mutation data/mut \
  --meta data/meta \
  --sitemask out/siteMask \
  --outDir out \
  --oJSON \
  --oGraphML \
  --nthread 3
```

The results will be generated in `example/out`. The pre-runned results would be found in `example/out_example`.

## Input Format

McAN can read variants from the Variant Call Format (VCF) or a mutation format. Here the mutation format is defined as follows:
```
<SampleName>\t<AccessionID>\t[<POS>(<VariantType>:<REF>-><ALT>)[;<POS>(<VariantType>:<REF>-><ALT>)]...]
```
SampleName can take missing value '*', in which case it will be read later from a metadata file. The mutation data must include the reference sample which contains no mutation (EPI_ISL_454904 in example below). Here's an example:
```
hCoV-19/England/ALDP-143B269/2021	EPI_ISL_1454756	1(Deletion:AATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGT->-);241(SNP:C->T);913(SNP:C->T);2062(SNP:C->T);3037(SNP:C->T);3267(SNP:C->T);5388(SNP:C->A);5986(SNP:C->T);6954(SNP:T->C);9746(SNP:C->T);11287(Deletion:GTCTGGTTTT->G);14408(SNP:C->T);14590(SNP:T->G);14676(SNP:C->T);15279(SNP:C->T);16176(SNP:T->C);16391(SNP:C->T);17615(SNP:A->G);19656(SNP:G->A);21764(Deletion:ATACATG->A);21990(Deletion:TTTA->T);23063(SNP:A->T);23271(SNP:C->A);23403(SNP:A->G);23604(SNP:C->A);23709(SNP:C->T);24506(SNP:T->G);24914(SNP:G->C);26461(SNP:C->T);26730(SNP:G->C);27972(SNP:C->T);28048(SNP:G->T);28111(SNP:A->G);28270(Deletion:TA->T);28280(SNP:G->C);28281(SNP:A->T);28282(SNP:T->A);28881(SNP:G->A);28882(SNP:G->A);28883(SNP:G->C);28973(SNP:A->G);28977(SNP:C->T);29109(SNP:C->A);29541(SNP:C->T);29640(SNP:C->T);29773(SNP:G->T);29836(Deletion:CCCATGTGATTTTAATAGCTTCTTAGGAGAATGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA->C)
hCoV-19/Denmark/DCGC-48559/2021	EPI_ISL_1124652	1(Deletion:AATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGT->-);158(SNP:A->G);241(SNP:C->T);913(SNP:C->T);3037(SNP:C->T);3267(SNP:C->T);5388(SNP:C->A);5986(SNP:C->T);6954(SNP:T->C);11287(Deletion:GTCTGGTTTT->G);12388(SNP:T->C);14408(SNP:C->T);14676(SNP:C->T);15279(SNP:C->T);16176(SNP:T->C);21764(Deletion:ATACATG->A);21990(Deletion:TTTA->T);23063(SNP:A->T);23271(SNP:C->A);23403(SNP:A->G);23604(SNP:C->A);23709(SNP:C->T);24506(SNP:T->G);24914(SNP:G->C);25448(SNP:A->G);26681(SNP:C->T);27972(SNP:C->T);28048(SNP:G->T);28111(SNP:A->G);28270(Deletion:TA->T);28280(SNP:G->C);28281(SNP:A->T);28282(SNP:T->A);28881(SNP:G->A);28882(SNP:G->A);28883(SNP:G->C);28977(SNP:C->T);29421(SNP:C->T);29422(SNP:G->T);29627(SNP:C->T);29836(Deletion:CCCATGTGATTTTAATAGCTTCTTAGGAGAATGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA->C)
hCoV-19/Germany/NW-RKI-I-063476/2021	EPI_ISL_1568686	1(Deletion:AAT->-);241(SNP:C->T);913(SNP:C->T);1937(SNP:T->C);1979(SNP:G->A);2110(SNP:C->T);3037(SNP:C->T);3267(SNP:C->T);4582(SNP:C->T);5388(SNP:C->A);5986(SNP:C->T);6954(SNP:T->C);7267(SNP:C->T);9856(SNP:G->T);11083(SNP:G->T);11287(Deletion:GTCTGGTTTT->G);14120(SNP:C->T);14408(SNP:C->T);14676(SNP:C->T);15279(SNP:C->T);16176(SNP:T->C);21766(Deletion:ACATGTC->A);21990(Deletion:TTTA->T);23063(SNP:A->T);23271(SNP:C->A);23403(SNP:A->G);23604(SNP:C->A);23690(SNP:A->G);23709(SNP:C->T);24506(SNP:T->G);24914(SNP:G->C);27972(SNP:C->T);28048(SNP:G->T);28095(SNP:A->T);28111(SNP:A->G);28270(Deletion:TA->T);28280(SNP:G->C);28281(SNP:A->T);28282(SNP:T->A);28881(SNP:G->A);28882(SNP:G->A);28883(SNP:G->C);28977(SNP:C->T);29732(Deletion:CCGAGGCCACGCGGAGTACGATCGAGTGTACAGTGAACAATGCTAGGGAGAGCTGCCTATATGGAAGAGCCCTAATGTGTAAAATTAATTTTAGTAGTGCTATCCCCATGTGATTTTAATAGCTTCTTAGGAGAATGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA->C)
hCoV-19/USA/AK-PHL8275/2021	EPI_ISL_1821911	1(Deletion:AATTAAAGGTTTATACCT->-);201(SNP:T->C);203(SNP:C->T);222(SNP:C->T);241(SNP:C->T);507(Deletion:ATGGTCATGTTATGGT->A);1738(SNP:G->T);3037(SNP:C->T);3140(SNP:C->T);9979(SNP:C->T);10029(SNP:C->T);10954(SNP:C->T);11117(SNP:A->G);12789(SNP:C->T);14408(SNP:C->T);17079(SNP:G->A);18255(SNP:G->T);19839(SNP:T->C);19974(SNP:A->G);21306(SNP:C->T);22995(SNP:C->A);23403(SNP:A->G);23604(SNP:C->A);23756(SNP:A->G);25553(SNP:C->T);28881(SNP:G->A);28882(SNP:G->A);28883(SNP:G->C);29197(SNP:C->T);29866(Deletion:ATGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA->A)
hCoV-19/USA/KS-KHEL-0753/2020	EPI_ISL_1412493	1(Deletion:AATTAAAGGTTTATACCTTCCCAG->-);219(SNP:G->T);241(SNP:C->T);3037(SNP:C->T);3122(SNP:G->T);4950(SNP:T->C);7108(SNP:C->T);7687(SNP:A->G);9409(SNP:A->G);11595(SNP:A->G);12970(SNP:C->T);14408(SNP:C->T);14712(SNP:G->T);16610(SNP:C->A);16733(SNP:C->T);18079(SNP:G->T);18555(SNP:C->T);19839(SNP:T->C);20069(SNP:G->T);23403(SNP:A->G);23426(SNP:G->T);23587(SNP:G->T);23756(SNP:A->G);24138(SNP:C->A);28881(SNP:G->A);28882(SNP:G->A);28883(SNP:G->C);29640(SNP:C->T);29870(Deletion:CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA->C)
hCoV-19/USA/TX-HMH-MCoV-18622/2020	EPI_ISL_784970	1(Deletion:AATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACC->-);241(SNP:C->T);1059(SNP:C->T);2388(SNP:C->T);2448(SNP:G->T);3037(SNP:C->T);4331(SNP:C->T);6040(SNP:C->T);6422(SNP:G->A);9112(SNP:C->M);10319(SNP:C->T);14408(SNP:C->T);15380(SNP:G->T);16926(SNP:T->C);18424(SNP:A->G);19974(SNP:A->G);21304(SNP:C->T);22225(SNP:G->T);23403(SNP:A->G);25563(SNP:G->T);25907(SNP:G->T);27964(SNP:C->T);28472(SNP:C->T);28869(SNP:C->T);29053(SNP:A->C);29095(SNP:C->T);29752(SNP:A->T);29854(Deletion:CTTCTTAGGAGAATGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA->C)
hCoV-19/USA/CA-CHLA-PLM57647613/2020	EPI_ISL_753540	241(SNP:C->T);1059(SNP:C->T);1927(SNP:T->C);3037(SNP:C->T);6196(SNP:C->T);10319(SNP:C->T);14408(SNP:C->T);15766(SNP:G->T);18424(SNP:A->G);18538(SNP:G->T);21304(SNP:C->T);22255(SNP:A->T);23120(SNP:G->T);23403(SNP:A->G);24146(SNP:C->G);25563(SNP:G->T);25907(SNP:G->T);25930(SNP:T->C);26951(SNP:G->T);27964(SNP:C->T);28472(SNP:C->T);28869(SNP:C->T);29439(SNP:A->T);29903(Insertion:A->AHMSHCAMTCSCYGWHMSHCVGSADBATCHCMSRSARSCVARAGTBACSSTSATCMR)
hCoV-19/Germany/SN-RKI-I-038205/2021	EPI_ISL_1284566	1(Deletion:AATTAAAGGTTTATACCTTCCCAGGTAACAA->-);241(SNP:C->T);815(SNP:C->T);913(SNP:C->T);1245(SNP:G->A);3037(SNP:C->T);3267(SNP:C->T);5388(SNP:C->A);5986(SNP:C->T);6954(SNP:T->C);11287(Deletion:GTCTGGTTTT->G);13446(SNP:C->T);14408(SNP:C->T);14676(SNP:C->T);15096(SNP:T->C);15279(SNP:C->T);16176(SNP:T->C);21764(Deletion:ATACATG->A);21990(Deletion:TTTA->T);23063(SNP:A->T);23271(SNP:C->A);23403(SNP:A->G);23604(SNP:C->A);23709(SNP:C->T);24506(SNP:T->G);24914(SNP:G->C);25785(SNP:G->T);27972(SNP:C->T);28048(SNP:G->T);28111(SNP:A->G);28270(Deletion:TA->T);28280(SNP:G->C);28281(SNP:A->T);28282(SNP:T->A);28881(SNP:G->A);28882(SNP:G->A);28883(SNP:G->C);28977(SNP:C->T);29764(SNP:G->A)
hCoV-19/Germany/un-RKI-I-145500/2021	EPI_ISL_2129747	1(Deletion:AATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCA->-);241(SNP:C->T);913(SNP:C->T);3037(SNP:C->T);3267(SNP:C->T);3302(SNP:G->A);3817(SNP:C->T);4002(SNP:C->T);5388(SNP:C->A);5944(SNP:C->T);5986(SNP:C->T);6954(SNP:T->C);11287(Deletion:GTCTGGTTTT->G);14408(SNP:C->T);14676(SNP:C->T);15096(SNP:T->C);15279(SNP:C->T);16176(SNP:T->C);21764(Deletion:ATACATG->A);21990(Deletion:TTTA->T);23063(SNP:A->T);23271(SNP:C->A);23403(SNP:A->G);23604(SNP:C->A);23709(SNP:C->T);24506(SNP:T->G);24914(SNP:G->C);25424(SNP:G->T);27389(SNP:C->T);27972(SNP:C->T);28048(SNP:G->T);28095(SNP:A->T);28111(SNP:A->G);28114(SNP:T->C);28270(Deletion:TA->T);28280(SNP:G->C);28281(SNP:A->T);28282(SNP:T->A);28881(SNP:G->A);28882(SNP:G->A);28883(SNP:G->C);28884(SNP:G->C);28977(SNP:C->T);29853(Deletion:GCTTCTTAGGAGAATGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA->G)
hCoV-19/Wuhan/HB-WH1-122/2020	EPI_ISL_454904
```
Metadata file is also required for building haplotype networks and assigning lineages, which contains 6 fields per record: SampleName, AccessionID, SamplingDate (yyyy-mm-dd or yyyy-mm), Country, State and City. All data lines are tab-delimited. Missing values are specified with a '*'. SampleName in metadata is only read and used when missing in the mutation data. The order of the samples in mutation data must match the order of the samples in metadata. Here's an example:
```
*	EPI_ISL_1454756	2021-03-13	United Kingdom	*	*
*	EPI_ISL_1124652	2021-02-15	Denmark	*	*
*	EPI_ISL_1568686	2021-03-22	Germany	*	*
*	EPI_ISL_1821911	2021-04-19	United States	*	*
*	EPI_ISL_1412493	2020-12-10	United States	*	*
*	EPI_ISL_784970	2020-11-20	United States	*	*
*	EPI_ISL_753540	2020-11-30	United States	*	*
*	EPI_ISL_1284566	2021-03-02	Germany	*	*
*	EPI_ISL_2129747	2021-04-30	Germany	*	*
*	EPI_ISL_454904	2020-03-02	China	*	*
```

## Options
```
Usage: McAN [option]...

Options:
  -f, --vcf <file>        input vcf file
  -u, --mutation <file>   input mutation file
  -m, --meta <file>       input meta file
  -s, --sitemask <file>   input sitemask file [optional]
  -x, --maxnsample <int>  maximum number of samples [default: Inf]
  -t, --nthread <int>     number of worker threads [default: 1]
  -o, --outDir <dir>      directory for output
  -J, --oJSON             output network with JSON format
  -T, --oTSV              output network with TSV format
  -G, --oGraphML          output network with GraphML format
  -M, --oMutation         convert vcf into mutation format and output

  -h, --help              display this help and exit
  -v, --version           output version information and exit


Usage: McAN siteMask [option]...

Options:
  -f, --vcf <file>       input vcf file
  -u, --mutation <file>  input mutation file
  -m, --meta <file>      input meta file
  -q, --minfreq <int>    filtered out variants with frequency < minfreq
  -o, --out <file>       output sitemask file
  -M, --oMutation        convert vcf into mutation format and output

  -h, --help             display this help and exit
  -v, --version          output version information and exit
```

## Contributors
- McAN contributors

## Citation
Please cite this for 2019 Novel Coronavirus Resource:
```
Zhao WM, Song SH, Chen ML, et al. The 2019 novel coronavirus resource. Yi Chuan. 2020;42(2):212–221. doi:10.16288/j.yczz.20-030 [PMID: 32102777]
```
or BibTeX Format:
```
@article{zhao20202019,
  title={The 2019 novel coronavirus resource.},
  author={Zhao, Wen-Ming and Song, Shu-Hui and Chen, Mei-Li and Zou, Dong and Ma, Li-Na and Ma, Ying-Ke and Li, Ru-Jiao and Hao, Li-Li and Li, Cui-Ping and Tian, Dong-Mei and others},
  journal={Yi chuan= Hereditas},
  volume={42},
  number={2},
  pages={212--221},
  year={2020}
}
```

## Contact
If you have any questions, feel free to contact us. lil@big.ac.cn

## License
MIT License
>  Copyright (c) 2022 McAN contributors
>
>  Permission is hereby granted, free of charge, to any person obtaining a copy
>  of this software and associated documentation files (the "Software"), to deal
>  in the Software without restriction, including without limitation the rights
>  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
>  copies of the Software, and to permit persons to whom the Software is
>  furnished to do so, subject to the following conditions:
>
>  The above copyright notice and this permission notice shall be included in
>  all copies or substantial portions of the Software.
>
>  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
>  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
>  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
>  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
>  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
>  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
>  THE SOFTWARE.
>
