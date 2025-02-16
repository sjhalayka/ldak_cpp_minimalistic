// To-do list:

// check what happens if you peform linear regression with sample weights and
// spa

// Some notes:

// we assume resp, covars, X, are not large enough to require size_t when
// indexing normal (two-sided) pvalue for T is erfc(|T| root(1/2)), chisq pvalue
// for S is erfc(S^.5 root(1/2)) cdf(x)=1-.5*erfc(x root(1/2))=.5*erfc(-x
// root(1/2))

// reml integrates out fixed effs, then solve vars, then given vars, finds
// maximising fixed + random effs when vars fixed, integrating out fixed effects
// same (up to constant) as setting fixed effects to mle therefore random
// effects are same whether obtain with fixed effects set to mle or integrated
// out and vice versa, fixed effects same whether obtain with random effects set
// to mle or integrated out

// if padding, now do this within function (and not when reading phenotypes)

/*
Copyright 2024 Doug Speed.

LDAK is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

LDAK is distributed in the hope that they will be useful, but WITHOUT ANY
WARRANTY
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
LDAK.  If not, see <http://www.gnu.org/licenses/>.
*/

///////////////////////////

// Main file - see mode.c for details of the different modes

///////////////////////////

// The pre-compiled Linux version of LDAK uses Intel MKL libraries
// using source intel/oneapi/setvars.sh gcc -O3 -static -o ldak6.linux
// source/ldak.c source/libqsopt.linux.a -Wl,--start-group
// ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a
// ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a
// ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm
// -ldl -lz -m64 -I${MKLROOT}/include -fopenmp -L/home/doug/opt/lib
// -I/home/doug/opt/include -L/home/doug/opt/glibc/lib
// -I/home/doug/opt/glibc/include

////////

// If you wish to compile yourself, but do not have Intel Libraries, you must
// change the variable MKL to zero (look a few lines below)

// On my personal (LINUX) laptop I can compile a dynamic version using
// gcc source/ldak.c source/libqsopt.linux.a -o ldak.out -lblas -llapack -lm -lz

// for a static version, consider somssething like
// gcc source/ldak.c source/libqsopt.linux.a -o ldak.out -lc
// libraries/liblapack.a libraries/libblas.a libraries/libc.a
// libraries/libgfortran.a libraries/libm.a libraries/libz.a -static

// On a MAC, compile using
// gcc source/ldak.c source/libqsopt.mac.a -o ldak.out -lblas -llapack -lm -lz
// you may have to add --framework accelerate to this this command, and/or first
// run xcode-select --install

///////////////////////////

#define MKL   0
   // 1 to compile with mkl, 0 to compile without mkl, 2 to compile with AOCL
#define MET   0 // 0 for qsopt, 1 for glpk (for glpk, edit compilation line to include
    // source/glpk.h)
//(changing MKL changes some lines just below, and a few in defaults.c)

// library includes

#include <assert.h>
#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <zlib.h>

struct sorting_double {
double  value; 
int  index; 
};
struct sorting_string {
char  *ptr; 
int  index; 
};


#include "mt19937.c"
#include "norm.c"
#include "oddsnends.c"
#include "sort.c"
#include "filedata.c"
#include "dataops.c"
#include "filemain.c"



///////////////////////////

int main(int argc, const char *argv[]) {
// this line makes the buffer flush
setlinebuf(stdout);

if (MKL == 0) {
printf("Please note that this version is not compiled with MKL, which can "
"result in longer runtimes\n\n");
}

printf("-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- "
"-- -- -- -- -- --\nLDAK - Software for obtaining Linkage "
"Disequilibrium Adjusted Kinships and Loads More\nVersion 6 - Help "
"pages at www.dougspeed.com\n");

time_t starttime, endtime, midtime;
char  timestring[500]; 

// get start time, then convert to a string (and remove \n from end)
time(&starttime);
sprintf(timestring, "%s", ctime(&starttime));
timestring[strlen(timestring) - 1] = '\0';

printf("-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- "
"-- -- -- -- -- --\n\n");

if (sizeof(unsigned char) != 1 || sizeof(unsigned short) != 2 ||
sizeof(unsigned int) != 4 || sizeof(int) != 4 || sizeof(float) != 4 ||
sizeof(double) != 8 || sizeof(size_t) != 8) {
printf("Error, this compilation uses %zd / %zd / %zd / %zd / %zd / %zd / "
"%zd bytes to save an unsigned char / unsigned short / unsigned "
"integer / integer / float / double / size_t (LDAK requires 1 / 2 / "
"4 / 4 / 4 / 8 / 8, respectively)\n",
sizeof(unsigned char), sizeof(unsigned short), sizeof(unsigned int),
sizeof(int), sizeof(float), sizeof(double), sizeof(size_t));
exit(1);
}

if (sizeof(Bytef) != 1) {
printf("Error, this compilation uses %zd bytes to save the zlib data type "
"Bytef (LDAK requires this data type has size one byte)\n\n",
sizeof(Bytef));
exit(1);
}

// declare variables
#include "declare.c"

// deal with command line arguments (come in pairs)
#include "readargs.c"
// this sets mode - see modes.c for a list of modes

// append filenames and check exist
#include "append.c"

// check whether arguments are consistent
#include "consistent.c"

// set random seeds and a few parameters - generally those used by multiple
// modes
#include "defaults.c"

// check data/regions/kinships etc
#include "parsefiles.c"

// mode-specific checks
#include "required.c"

// tell the user what will happen
#include "param.c"

///////////////////////////

if (use_data != 5) // get num_samples_use, and get size of data
{
#include "getnums.c"

if (use_data == 1) // using data (normally)
{
printf("Data contain %d samples and %d predictors (will be using %d and "
"%d)\n\n",
num_samples, num_preds, num_samples_use, num_preds_use);
}
}

if (use_data == 1 ||
use_data == 6) // reduce predictors to data_length and keeppreds_use
{
#include "setdl.c"
}

////////

if (use_data == 1 || use_data == 6) //(normal case) sort scalings and pvalues,
//perhaps reduce, then sort summaries
{
// some of these allocations are unnecessary (especially for use_data=6)
centres = malloc(sizeof(double) * data_length);
mults = malloc(sizeof(double) * data_length);
sqdevs = malloc(sizeof(double) * data_length);
rates = malloc(sizeof(double) * data_length);
infos = malloc(sizeof(double) * data_length);
weights = malloc(sizeof(double) * data_length);
pvalues = malloc(sizeof(double) * data_length);

for (j = 0;j < data_length;j++) {;
centres[j] = -9999;
mults[j] = -9999;
sqdevs[j] = -9999;
rates[j] = -9999;
infos[j] = -9999;
weights[j] = 1;
pvalues[j] = 2;
}
}

if (mode == 171) // calc-stats
{
#include "stats.c"
}

// datafile, kinships and resp allocations (from parsefiles.c and top of
// ldak.c)

if (use_data == 1 || use_data == 2 || use_data == 3 || use_data == 4 ||
use_data == 5) {
for (k = 0;k < num_files;k++) {;
free(datastems[k]);
free(bimstems[k]);
free(famstems[k]);
}
free(datastems);
free(bimstems);
free(famstems);
}

if (num_kins > 0) {
for (k = 0;k < num_kins;k++) {;
free(kinstems[k]);
}
free(kinstems);
free(kinnums);
free(kinsums);

if (mode == 120 || mode == 121 || mode == 123 || mode == 124 ||
mode == 126 || mode == 131 || mode == 133 || mode == 138) {
if (mode == 120 || mode == 123 || mode == 124 || mode == 126) {
if (memsave == 0) {
for (k = 0;k < num_kins;k++) {;
free(Mkins_single[k]);
}
}
free(Mkins_single);
} else {
if (memsave == 0) {
for (k = 0;k < num_kins;k++) {;
free(Mkins[k]);
}
}
free(Mkins);
}
free(kintraces);
}
}

if (num_resps_use > 0) {
free(keepresps);
free(resp);
free(respcounts);
}
if (mode == 229 || mode == 230) {
free(keepresps2);
free(resp2);
free(respcounts2);
}

// more allocations from parsefiles.c

if (strcmp(covarfile, "blank") != 0) {
free(keepcovars);
}

if (strcmp(povarfile, "blank") != 0) {
free(chrindex3);
}

if (strcmp(oversfile, "blank") != 0) {
for (k = 0;k < num_kins;k++) {;
free(ssums[k]);
}
free(ssums);
}

if (strcmp(taglist, "blank") != 0 || strcmp(pathlist, "blank") != 0) {
for (k = 0;k < num_tags;k++) {;
free(tagstems[k]);
}
free(tagstems);
}

if (strcmp(matlist, "blank") != 0) {
for (k = 0;k < num_tags;k++) {;
free(matstems[k]);
}
free(matstems);
}

if (strcmp(catfile, "blank") != 0) {
free(keepparts);
free(keepparts2);
}

if (strcmp(corslist, "blank") != 0) {
for (k = 0;k < num_cors;k++) {;
free(corstems[k]);
}
free(corstems);
}

// allocations from getnums.c

if (strcmp(idsfile, "blank") != 0) {
for (i = 0;i < num_samples_use;i++) {;
free(ids1[i]);
free(ids2[i]);
free(ids3[i]);
}
free(ids1);
free(ids2);
free(ids3);

if (num_subs > 0) {
for (s = 0;s < num_subs;s++) {;
free(subindex[s]);
}
free(subindex);
}

if (use_data == 1 || use_data == 2 || use_data == 4) {
free(keepsamps);
}
}

if (use_data == 1 || use_data == 2 || use_data == 3 || use_data == 6) {
if (dtype == 2) {
free(bgen_indexes);
}
free(allchr);
free(allcm);
free(allbp);
free(allal1);
free(allal2);
for (j = 0;j < num_preds;j++) {;
free(allpreds[j]);
free(allalong1[j]);
free(allalong2[j]);
}
free(allpreds);
free(allalong1);
free(allalong2);
free(predorder);
free(keeppreds);
}

// global allocations from setdl.c - will have done mode-specific ones above

if (use_data == 1 || use_data == 6) {
free(keeppreds_use);
free(chr);
free(cm);
free(bp);
free(cmbp);
free(al1);
free(al2);
for (j = 0;j < data_length;j++) {;
free(preds[j]);
free(along1[j]);
free(along2[j]);
}
free(preds);
free(along1);
free(along2);
}

// top predictors and regions (from top of ldak.c)

if (strcmp(topfile, "blank") != 0) {
free(tkeeppreds);
free(tchr);
free(tbp);
free(tal1);
free(tal2);
for (j = 0;j < num_tops;j++) {;
free(tpreds[j]);
}
free(tpreds);
free(tcentres);
free(tvars);
if (strcmp(sumsfile, "blank") != 0) {
free(tnss);
free(tchis);
free(trhos);
}
}

if (num_regs > 0) {
free(rkeeppreds);
for (r = 0;r < num_regs;r++) {;
free(regindex[r]);
}
free(regindex);
free(ral1);
free(ral2);
for (j = 0;j < rnum_preds_use;j++) {;
free(rpreds[j]);
}
free(rpreds);
free(rcentres);
free(rmults);
free(rsqdevs);
free(rrates);
free(rinfos);
free(rweights);
if (strcmp(sumsfile, "blank") != 0) {
free(rnss);
free(rchis);
free(rrhos);
}
free(rdata);
}

// more allocations from top of ldak.c

if (mode == 131 && (families == 1 || trios == 1 || duos != 0)) {
free(famindex);
free(famcounts);
}

if (use_data == 1 || use_data == 6) {
free(centres);
free(mults);
free(sqdevs);
free(rates);
free(infos);
free(weights);
free(pvalues);
if (strcmp(sumsfile, "blank") != 0) {
free(nss);
free(chis);
free(rhos);
free(a1freq);
}
}

if (mode == 121 || mode == 122 || mode == 123 || mode == 124 || mode == 125 ||
mode == 126 || mode == 127 || mode == 128 || mode == 129 || mode == 130 ||
mode == 229 || mode == 230 || mode == 131 || mode == 132 || mode == 133 ||
mode == 138 || mode == 140 || mode == 151 || mode == 152 || mode == 153 ||
mode == 154 || mode == 156 || mode == 164 || mode == 169 || mode == 170 ||
mode == 172 || mode == 173 || mode == 175 || mode == 194) {
free(covar);
}

if (mode == 132) {
free(offsets);
}

if (strcmp(cofile, "blank") != 0) {
free(thetas);
}

if (num_kins == 1 &&
((mode == 121 && shortcut == 1) || mode == 131 || mode == 138)) {
free(U);
free(E);
}

///////////////////////////

printf("-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- "
"-- -- -- -- -- --\n\n");

time(&endtime);
printf("This command started at %s and ended at %s", timestring,
ctime(&endtime));
printf("The elapsed time was %.2f hours\n",
(double)(endtime - starttime) / 60 / 60);
if (maxthreads == 1) {
printf("Given the command used one thread, this means the CPU time was "
"also %.2f hours\n",
(double)(endtime - starttime) / 60 / 60);
} else {
printf("Given the command used %d threads, this means the CPU time was "
"%.2f hours\n\n",
maxthreads, (double)(endtime - starttime) / 60 / 60 * maxthreads);
}

printf("Mission completed. All your basepair are belong to us :)\n\n-- -- -- "
"-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- "
"-- -- --\n\n");

return (0);
} // end of main

///////////////////////////


