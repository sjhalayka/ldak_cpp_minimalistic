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

// Declare variables 

/////////////////////////// 

// parameters set directly by user 

int  dougvar = 0;
double  dougvar2 = -9999;

// mode indicates which main argument is called (see mode.c for list) 
int  mode = -9999;

// random seed 
int  seed = -9999;

// data input 

char  folder2[500] = "blank";
char  folder[500] = "";
char  outfile2[500] = "blank";
char  outfile[500] = "";

char  workdir2[500] = "blank";
char  workdir[500] = "";

char  udatafile2[500] = "blank";
char  udatafile[500] = "";
char  ubimfile2[500] = "blank";
char  ubimfile[500] = "";
char  ufamfile2[500] = "blank";
char  ufamfile[500] = "";
char  unamefile2[500] = "blank";
char  unamefile[500] = "";
char  datalist2[500] = "blank";
char  datalist[500] = "";

int  oxchr = -9999;
int  genskip = -9999;
int  genheaders = -9999;
int  genprobs = -9999;
int  nonsnp = -9999;

int  bedzeros[256];
int  bedones[256];
int  bedtwos[256];
double  missingvalue = -9999;

// data filtering 

char  bsampfile2[500] = "blank";
char  bsampfile[500] = "";
char  csampfile2[500] = "blank";
char  csampfile[500] = "";
int  num_subs = -9999;
char  subpref2[500] = "blank";
char  subpref[500];

char  bpredfile2[500] = "blank";
char  bpredfile[500] = "";
char  cpredfile2[500] = "blank";
char  cpredfile[500] = "";
int  onechr = -9999;
char  onesnp[500] = "blank";

double  minmaf = -9999;
double  maxmaf = -9999;
double  minvar = -9999;
double  minobs = -9999;
double  mininfo = -9999;

// data scaling (and pvalues) and coding 

char  centresfile2[500] = "blank";
char  centresfile[500] = "";
char  weightsfile2[500] = "blank";
char  weightsfile[500] = "";
int  ignoreweights = -9999;
double  power = -9999;
int  hwestand = -9999;
char  pvafile2[500] = "blank";
char  pvafile[500] = "";
char  impfile2[500] = "blank";
char  impfile[500] = "";

int  encoding = -9999;
double  threshold = -9999;
double  minprob = -9999;

// kinships, regions, responses, summaries and fixed 

char  kinname2[500] = "blank";
char  kinname[500] = "";
char  kinlist2[500] = "blank";
char  kinlist[500] = "";
int  kindetails = -9999;

int  num_regs = -9999;
char  regpref2[500] = "blank";
char  regpref[500];
double  rprune = -9999;

char  respfile2[500] = "blank";
char  respfile[500] = "";
int  mpheno = -9999;
int  mpheno2 = -9999;
int  pad = -9999;

char  sumsfile2[500] = "blank";
char  sumsfile[500] = "";
char  sums2file2[500] = "blank";
char  sums2file[500] = "";

int  fixn = -9999;
int  fixn2 = -9999;
int  amb = -9999;
double  scaling = -9999;
double  scaling2 = -9999;

double  prev = -9999;
double  prev2 = -9999;
double  ascer = -9999;

char  covarfile2[500] = "blank";
char  covarfile[500] = "";
char  covarnums2[2000] = "blank";
char  covarnums[2000] = "blank";
char  covarnames2[2000] = "blank";
char  covarnames[2000] = "blank";
char  envfile2[500] = "blank";
char  envfile[500] = "";
char  topfile2[500] = "blank";
char  topfile[500] = "";
char  factorfile2[500] = "blank";
char  factorfile[500] = "";
char  povarfile2[500] = "blank";
char  povarfile[500] = "";

char  offsetfile2[500] = "blank";
char  offsetfile[500] = "";

// calculating weights, thinning and finding/removing tags 

int  nothin = -9999;
double  wprune = -9999;
double  window_kb = -9999;
double  window_kb2;
int  window_length = -9999;
int  window_length2;
double  window_cm = -9999;

double  section_kb = -9999;
int  section_length = -9999;
double  section_cm = -9999;
double  buffer_kb = -9999;
int  buffer_length = -9999;
double  buffer_cm = -9999;

int  section = -9999;
int  section_start = -9999;

int  lddecay = -9999;
double  halflife = -9999;

int  fudge = -9999;
int  simplex = -9999;
double  maxtime = -9999;
int  spread = -9999;

char  targetfile2[500] = "blank";
char  targetfile[500];

// calculating and manipulating kinships (partitions also used for 
// gene/chunk-based reml, sumher, gre) 

int  part_length = -9999;
int  bychr = -9999;
int  num_parts = -9999;
char  partpref2[500] = "blank";
char  partpref[500];
int  checkpart = -9999;

int  partition = -9999;
int  kingz = -9999;
int  kinraw = -9999;
int  single = -9999;

double  dosage = -9999;
char  malesfile2[500] = "blank";
char  malesfile[500];
int  onlydets = -9999;
char  invsfile2[500] = "blank";
char  invsfile[500];
int  david = -9999;

double  maxrel = -9999;
double  minrel = -9999;
int  kinstand = -9999;

int  partial = -9999;

// reml, blup and he/pcgc (shortcut used indirectly by other functions) 

int  diagonal = -9999;
char  hersfile2[500] = "blank";
char  hersfile[500] = "";
int  hestart = -9999;
int  shortcut = -9999;
int  discenv = -9999;
char  oversfile2[500] = "blank";
char  oversfile[500] = "";

char  remlfile2[500] = "blank";
char  remlfile[500] = "";

int  adjusted = -9999;
double  trun = -9999;

int  num_vects = -9999;
int  ldlt = -9999;

char  relfile[500];
char  relfile2[500] = "blank";
int  cordups = -9999;

// association analysis (genes also used for condensing data) 

int  kvikstep = -9999;
int  faststep = -9999;
int  gctastep = -9999;
char  prsfile2[500] = "blank";
char  prsfile[500] = "";
int  verbose = -9999;

int  spatest = -9999;
int  num_knots = -9999;
int  num_bins = -9999;
int  spaside = 1;
double  spathresh = -9999;
double  spamax = -9999;

int  families = -9999;
int  trios = -9999;
int  duos = -9999;

int  adjpreds = -9999;
char  sampwfile2[500] = "blank";
char  sampwfile[500] = "";
int  sandwich = -9999;
int  exact = -9999;
int  scoretest = -9999;
char  transfile2[500] = "blank";
char  transfile[500] = "";

char  genefile2[500] = "blank";
char  genefile[500] = "";
double  chunks = -9999;
int  chunksbp = -9999;
int  gene_buffer = -9999;
int  up_buffer = -9999;
int  down_buffer = -9999;
double  minweight = -9999;
int  overlap = -9999;

int  gene_perms = -9999;
int  saveall = -9999;
double  gprune = -9999;
double  limit = -9999;

int  magma = -9999;
double  cut1 = -9999;
double  cut2 = -9999;
double  gamp = -9999;
double  gam1 = -9999;
double  gam2 = -9999;

// sumher (annotations also used by fast he and pcgc 

int  num_anns = -9999;
char  annpref2[500] = "blank";
char  annpref[500];
char  genpreds2[500] = "blank";
char  genpreds[500];
char  labfile2[500] = "blank";
char  labfile[500];
int  backpart = -9999;
int  allone = -9999;
int  reduce = -9999;

char  printfile2[500] = "blank";
char  printfile[500] = "";
char  herfile2[500] = "blank";
char  herfile[500] = "";
int  unbias = -9999;
int  savemat = -9999;
int  cover = -9999;
int  fourdp = -9999;

char  taglist2[500] = "blank";
char  taglist[500] = "";
char  matlist2[500] = "blank";
char  matlist[500] = "";
char  pathlist2[500] = "blank";
char  pathlist[500] = "";
int  checkdups = -9999;

char  tagfile2[500] = "blank";
char  tagfile[500] = "";
char  pathfile2[500] = "blank";
char  pathfile[500] = "";
char  altfile2[500] = "blank";
char  altfile[500] = "";
char  cvsfile2[500] = "blank";
char  cvsfile[500] = "";
char  catfile2[500] = "blank";
char  catfile[500] = "";
char  taufile2[500] = "blank";
char  taufile[500] = "";

int  checksums = -9999;
int  gcon = -9999;
int  cept = -9999;

int  ldsc = -9999;
int  chisol = -9999;
int  tagone = -9999;

int  divide = -9999;
int  uptaus = -9999;
char  powfile2[500] = "blank";
char  powfile[500] = "";

int  plet = -9999;

char  matfile2[500] = "blank";
char  matfile[500] = "";

char  expfile2[500] = "blank";
char  expfile[500] = "";
double  cvar = -9999;

// individual-level data prediction, then megaprs 

char  indhers2[500] = "blank";
char  indhers[500] = "";
double  herscale = -9999;

int  loco = -9999;
int  dichot = -9999;
int  multi = -9999;
int  fast = -9999;
int  fprs = -9999;
int  fastgwa = -9999;
char  fastfile2[500] = "blank";
char  fastfile[500] = "";

int  skipcv = -9999;
double  cvprop = -9999;
char  bvsfile2[500] = "blank";
char  bvsfile[500] = "";

int  ldpred = -9999;
int  pointmass = -9999;

int  ndivs = -9999;
int  nmcmc = -9999;
double  maxher = -9999;

int  checkped = -9999;
int  nped = -9999;
double  maithresh = -9999;
double  maisig = -9999;
int  ncal = -9999;
int  ncomp = -9999;
double  cthresh = -9999;
int  nscan = -9999;
int  revher = -9999;

char  blockfile2[500] = "blank";
char  blockfile[500] = "";

char  corslist2[500] = "blank";
char  corslist[500] = "";

double  subprop = -9999;

char  corname2[500] = "blank";
char  corname[500] = "";
char  pseudostem2[500] = "blank";
char  pseudostem[500] = "";

int  ptype = -9999;
char  bestfile2[500] = "blank";
char  bestfile[500] = "";

int  checkld = -9999;
char  ldfile2[500] = "blank";
char  ldfile[500] = "";

char  fracfile2[500] = "blank";
char  fracfile[500] = "";

int  checkfreq = -9999;
int  prsvar = -9999;
double  jackprop = -9999;

// pca, decompose, adjust-grm and others 

int  axes = -9999;
char  pcastem2[500] = "blank";
char  pcastem[500] = "";

int  eigenraw = -9999;
char  eigenfile2[500] = "blank";
char  eigenfile[500] = "";

int  noise = -9999;

// stats, scores, making phenotypes, snps, inflation, jackknifing, folds, find 
// gaussian 

char  scorefile2[500] = "blank";
char  scorefile[500] = "";
char  cofile2[500] = "blank";
char  cofile[500] = "";
int  savecounts = -9999;
char  finalfile2[500] = "blank";
char  finalfile[500] = "";

int  num_phenos = -9999;
int  num_causals = -9999;
double  her = -9999;
double  cher = -9999;
double  bivar = -9999;
double  bivar2 = -9999;
double  bivar3 = -9999;
char  probsfile2[500] = "blank";
char  probsfile[500];
char  causalsfile2[500] = "blank";
char  causalsfile[500];
char  effectsfile2[500] = "blank";
char  effectsfile[500];

int  num_inds = -9999;
int  num_snps = -9999;
double  maf1 = -9999;
double  maf2 = -9999;
int  nchrom = -9999;
int  famsize = -9999;
double  closeness = -9999;
int  quads = -9999;
int  pops = -9999;
double  missrate = -9999;

char  predlista[500];
char  predlista2[500] = "blank";
char  predlistb[500];
char  predlistb2[500] = "blank";
int  savepairs = -9999;

char  jackfile[500];
char  jackfile2[500] = "blank";
char  proffile[500];
char  proffile2[500] = "blank";
int  auc = -9999;

int  num_folds = -9999;

char  likefile[500];
char  likefile2[500] = "blank";
int  num_means = -9999;
int  num_sds = -9999;
double  minmean = -9999;
double  maxmean = -9999;
double  maxsd = -9999;
int  omitone = -9999;

// making and condensing data 

int  exsame = -9999;
int  exlong = -9999;
int  speedlong = -9999;

int  useminor = -9999;

// gre 

int  sinv = -9999;
char  greout[500];
char  greout2[500] = "blank";

// common options 

int  checkroot = -9999;
double  mincor = -9999;
double  maxcor = -9999;
double  cutoff = -9999;
int  constrain = -9999;
int  num_blocks = -9999;
int  permute = -9999;
double  shrink = -9999;
double  strip = -9999;

int  bitsize = -9999;
int  bitsize2;
double  tol = -9999;
int  maxiter = -9999;
int  memsave = -9999;

int  manysamples = -9999;
int  manypreds = -9999;

// threading option 

int  maxthreads = 1;

/////////////////////////// 

// variables which are used for a specific purpose (and set fairly early on) 

// storing data 

int dtype = -9999, binary = -9999, famhead = -9999, extract = 0, use_data,
    num_files = 0;
char  **datastems;
char  **bimstems;
char  **famstems;
char  datafile[500];
char  bimfile[500];
char  famfile[500];

int num_samples = -9999, num_samples_use = -9999, *keepsamps, *keepsamps_males,
    **subindex;
int num_preds = -9999, num_preds_use = -9999, *keeppreds, data_length = -9999,
    *keeppreds_use;

int bgen_preds = -9999, bgen_samples = -9999, bgen_comp = -9999,
    bgen_layout = -9999, bgen_ids = -9999;
size_t  *bgen_indexes;
int  maxpreds = 500000;
// the apriori max number of predictors in gen files 

char  **ids1;
char  **ids2;
char  **ids3;
char  **ids4;
char  **allids1;
char  **allids2;
char  **allids3;
int *chr, *allchr, num_chr, num_chr2, num_chr3, *chrindex, *chrindex2,
    *chrindex3;
double  *cm;
double  *bp;
double  *cmbp;
double  *allcm;
double  *allbp;
double  *chrprops;
char  **preds;
char  **allpreds;
char  **along1;
char  **along2;
char  **allalong1;
char  **allalong2;
char  *al1;
char  *al2;
char  *allal1;
char  *allal2;

unsigned  char **data_char;
float  *data_single;
float  *data_single2;
double  *data;
double  *data2;
double  *data3;
double  **bytelookup;

float  *speedstarts;
float  *speedscales;
float  **ps;
float  **ps2;
double  *centres;
double  *mults;
double  *sqdevs;
double  *rates;
double  *infos;
double  *weights;
double  *pvalues;

int  num_pows;
double  *powers;
double  *powers2;

// kinships, regions, responses, summaries and fixed 

int  num_kins;
int  *kinnums;
float  *kins_single;
float  *kins_single2;
float  **Mkins_single;
float  *inv_single;
double  *kins;
double  *kins2;
double  **Mkins;
double  *kintraces;
double  *kinsums;
char  **kinstems;

int  *kindex;
double  *kcentres;
double  *kmults;
double  *kweights;
double  *kexps;
char  **kpreds;
char  *kal1;
char  *kal2;

int  rnum_preds_use;
int  *rkeeppreds;
int  **regindex;
int  rdata_length;
double  *rdata;
double  *rcentres;
double  *rmults;
double  *rsqdevs;
double  *rrates;
double  *rinfos;
double  *rweights;
char  **rpreds;
char  *ral1;
char  *ral2;

int num_resps, num_resps_use, *keepresps, *keepresps2, *respcounts,
    *respcounts2;
double  *resp;
double  *resp2;

int  sformat;
int  sformat2;
int  gotfreq;
double  *nss;
double  *nss2;
double  *nss3;
double  *rnss;
double  *tnss;
double  *snss;
double  *snss2;
double  *chis;
double  *chis2;
double  *chis3;
double  *rchis;
double  *tchis;
double  *schis;
double  *schis2;
double  *rhos;
double  *rhos2;
double  *rhos3;
double  *rrhos;
double  *trhos;
double  *srhos;
double  *srhos2;
double  *a1freq;

int  num_quants;
int  num_cats;
int  num_covars;
int  *keepcovars;
int  num_envs;
int  num_fixed;
int  num_prs;
double  *covar;
double  *thetas;
double  *thetasds;
double  *thetapvas;
double  covher;
double  topher;

int  num_tops;
int  *tkeeppreds;
int  *tchr;
double  *tbp;
double  *tcentres;
double  *tvars;
char  **tpreds;
char  *tal1;
char  *tal2;

double  *offsets;

// calculating weights, thinnings and finding/removing tags 

int  num_sections;
double  decay;
int  *sstarts;
int  *sends;
int  *sstarts2;
int  *sends2;
float  *cors_single;
double  *cors;
double  *cors2;
double  *datasqs;
double  *tally1;
double  *tally2;
double  *tally3;
double  *tally4;
int  *replace;

int  num_seek;
int  *sindex;
int  *stypes;

// calculating and manipulating kinships 

int  *pstarts;
int  *pends;
int  *highs;
int  *losts;
double  *maxes;

// reml, blup and he/pcgc 

char  blupfile2[500] = "blank";
char  blupfile[500] = "";
char  regfile2[500] = "blank";
char  regfile[500] = "";

int  **nums;
double **blupcentres, **blupfactors, *bluprands, *blupvalues, *blupmids,
    *blupprojs;
double  **guesses;
double  ***guesses2;

int  *bstarts;
int  *bends;
double  *R;
double  *RTdata;
double  *RTdata2;
double  *MkinsD;
double  *MkinsD2;
double  *MkinsR;
double  *MkinsR2;
double  *KKtraces;
double  *KKtraces2;
double  *KKtraces3;
double  *KYtraces;
double  *KYtraces2;

float  *bluprand_single;
float  *guess_single;

// association analysis 

int  num_fams;
int  *famindex;
int  *famcounts;

double *Pincs, *Pscales, *Pscales2, *Pscales3, *Pvarphens, *Plambdas, *Proots,
    *Pspamax;

int  *tindex;
int  *spastatus;
double  *knots;
double  *bins;
double  **CGF0;
double  **CGF1;
double  **CGF2;
double  **CGF3;
double cuts[6] = {.1, .01, .001, .0001, .00001, 5e-8};
double  *sweights;
double  DTD[9];
double  DTD2[3];
double  DTD3[4];
double  DTS[3];
double  DTS2[3];
double  invDTD[4];
double  *nullprobs;
double  *nullprobs2;
double  *nullweights;
double  *nullweights2;

int  kvikparity;
double  *prs;
double  *covprs;

int  num_genes;
int  genemax;
int  *gchr;
int  *gstrand;
int  *gstarts;
int  *gends;
int  *gparts;
double  *gbp1;
double  *gbp2;
double  *gpvas;
char  **gnames;
int  nclump;
double  *permlike;

// tagging, sumher and bayes factors 

int  **pindexes;
int  addpart;
int  addpart2;
int  addgenic;
int  *windex;
int  *vindex;
double  **pweights;

int  num_tags;
int  num_anals;
int  num_reds;
int  parttype = -9999;
char  **tagstems;
char  **matstems;
int  *pcounts;
int  *keepparts;
int  *keepparts2;

double  *exps;
double  *exps2;
char  **catlabels;

int  num_parts_use;
int  ncv;
int  *cvindex;
double  *stags;
double  **svars;
double  **svars2;
double  **ssums;
double  **ssums2;
double  **ssums_use;
double  *cvexps;
char  **spreds;
char  *sal1;
char  *sal2;

int  topm;
int  tops;
double  *influs;
double  *powlikes;
double  *mvals;
double  *svals;
double  **perfs;

double  *pmeans;
double  *pvars;
double  *bfs;
double  *podds;

// individual-level data prediction, then megaprs 

int  num_train;
int  num_test;
int  *keeptrain;
int  *keeptest;
int  num_try;
int  num_hers1;
int  num_hers2;
int  num_small;
int  num_full;
int  useped;
int  highstruct;
int  usecal;
int  usecomp;

int  *pedtops;
double  *pedscales;
double  *pedhers;
double  *pedprs;
double  *pedgammas;
double  *pedsds;
double  *pedmse;

int  *eindex;
int  *eindex2;
double  *ecentres;
double  *emults;
double  *esqdevs;
double  *erates;
double  *einfos;

int  *cindex;
int  *cindex2;
double  *cdata;
double  *ccentres;
double  *cmults;
double  *csqdevs;
double  *crates;
double  *cinfos;

int  *dindex;
int  *dindex2;
double  *ddata;
double  *dcentres;
double  *dmults;
double  *dsqdevs;
double  *drates;
double  *dinfos;

int  maxpairs;
int  num_rels;
int  *firsts;
int  *seconds;
double  *kinships;
double  *cX;
double  *cR;

double  *lambdas;
double  *lambdas2;
double  *lambdas3;
double  *lambdas4;
double  *cgammas;
double  *csds;
double  *ceffs;
double *effs, *effs2, *effs3, *probs, *probs2, *residuals, *residuals2,
    *residuals3, *changes, *ess, *ess2;

int  *bitrun;
int  *bitdo;
int  *bitactive;
int  *bitdet1;
int  *bitdet2;
double  *bitdiffs;
double  *bitpens;

int  *Mtops;
int  *Mbests;
int  *Mincs;
double  *Mscales;
double  *Mmses;
double  *Mmses2;
double  *Mneffs;
double *Gmat, *Emat, *Umat, *Umat2, *Wmat, *Wmat2, *Wmat3, *Wmat4, *Wmat5,
    *Vmat, *Vmat2, *Ymat;

// double *cY, *cX, *cR, *cP, *cHP, *cKP, *cVP, *ctemps, *calphas, *cbetas 
// double *results1a, *results1b, *results2a, *results2b 

size_t  *blockindexes;
int  *blockstarts;
int  *blockends;
int  *blocklengths;
int  **blockuse;
double  *rjksums;
double  *rjkaves;
double  *rjktemp;
double  *randnorms;
double  *randnorms2;
double  *datarands;
double  *datarands2;
double  *subrhos;
double  *restrhos;

int  num_cors;
char  **corstems;

int  *trytypes;
double *tryhers, *trylams, *tryscales, *tryps, *tryp2s, *tryp3s, *tryp4s,
    *tryf2s, *tryvars;

int  *highlds;
double  *loads;
double  *loads2;

// pca, decompose, adjust-grm and others 

float *Z_single, *ZTZ_single, *kinZ_single, *ZTkinZ_single, *W_single,
    *WZTkinZ_single;

// stats, scores, making phenotypes, snps, inflation, folds, jackknifing and 
// cors 

int  *presents;
int  *hets;
int  num_scores = 0;
double  *indsds;
double  *predmeans;
double  *predvars;
double  *predcors;
double  *predsds;
double  **unders;
double  **phens;
double  *vsds;
double  *veffs;
double  *vadds;
double  **liabs;

int  numa;
int  *keepa;
int  numb;
int  *keepb;
double *acentres, *amults, *bcentres, *bmults, *asqdevs, *bsqdevs, *arates,
    *brates, *ainfos, *binfos;
char  **apreds;
char  **bpreds;

int  *sZ1;
int  *sZ2;
double  *sX;
double  *sXTX;
double  *sXTX2;
double  *sW;
double  *sT;

// making and condensing data 

int  *Xcurrent;
int  *Xnall;
int  *Xnuse;
int  **Xks;
int  **Xks2;
int  *XNall;
int  *XNuse;
int  **Xkp;
int  **Xkp2;
char  ***Xids1;
char  ***Xids2;
char  ***Xids3;
int  **Xchr;
double  **Xcm;
double  **Xbp;
char  ***Xpreds;
char  ***Xalong1;
char  ***Xalong2;
char  **Xal1;
char  **Xal2;
FILE  **Xinput;
gzFile  *Xdatainputgz;

int  passqc;
int  csamps;
int  qmerge;
int  erra;
int  errb;
int  errc;
int  errd;
char  **pid;
char  **mid;
char  **bid;
char  **schar;
char  **pchar;
unsigned  char onechar;
unsigned  short oneshort;
float  minfloat;
float  maxfloat;

int  rowlength;
int  c0;
int  c1;
int  c2;
unsigned  char startchars[3];
unsigned  *rowchars;
unsigned  *rowchars2;

// speed tests 

float  *Rvsing;
float  *Rvsing2;
float  *Rvsing3;
float  *kins_packed;
double  *Rv;
double  *Rv2;
size_t  scount2;

// common options 

int bit, bittotal, bittotal2, bitstart, bitstart2, bitend, bitend2, bitlength,
    bitlength2, bitmax, step;

double  *Y;
double  *Yadj;
double  *Yadj2;
double  *YTdata;
double  *YTdata2;
double  YTCY;
double  *Z;
double  *Z2;
double  *Z3;
double  *ZTZ;
double  *ZTZ2;
double  *ZTZ3;
double  detZTZ;
double  *ZTdata;

double  *U;
double  *U2;
double  *E;
double  *E2;
double  *UTY;
double  *UTZ;
double  *UTdata;

int  Xtotal;
int  *Xstarts;
int  *Xends;
double  *X;
double  *XTCX;
double  *XTCX2;
double  *XTCX3;
double  *XTCX4;
double  *Xbeta;
double  *Xsums;
double  *Xnss;
double  *Xrhos;
double  *Xsqs;

// threading options 

int  thread;
int  threadstart;
int  threadend;
int  *Mcurrent;
float  **Mdatatemp;
FILE  **Minput;

/////////////////////////// 

// generic working variables 

size_t  scount;
size_t  smax;
int i, i2, i3, j, j2, j3, k, k2, g, m, m2, p, p2, q, q2, q3, r, s, count,
    count2, count3, count4;
int current, head, found, total, total2, total3, token, indcount, ecount,
    wcount, xcount, *ycounts;
int shuffle, start, end, best, worst, mark, mark2, mark3, mark4, gen, gen2,
    *gens, flag, cflag, eflag, hflag, *order, *order2, cols[6];
double sum, sum2, sum3, sumsq, sumsq2, sumsq3, sumsq4, mean, mean2, mean3, var,
    var2, var3, value, value2, value3, value4, value5;
double  last;
double  min;
double  max;
double  med;
double  maf;
double  varphen;
double  gif;
double  postmean;
double  unifrand;
double *hers, *hersold, *hersds, *shares, *sharesds, *cohers, *cohers2, mat[4],
    mat2[2];
double  likenull;
double  like;
double  lrtstat;
double  lrtpva;
double  *likes;
double  *likesold;
double  neff;
double  neff2;
double  neff3;
double  *pens;
double  factor;
double  weightsum;
double  minpvalue;
double  minpvalue2;
double  *jacks;
double  *gaussian;
double  *polates;
double  *polates2;
double  *vstarts;
double  *stats;
double  *stats2;
double  *stats3;
double  *stats4;
double  **effects;
double  *varexp;

int  one = 1;
int  two = 2;
int  three = 3;
int  info;
int  lwork;
int  *iwork;
int  *ipiv;
int  *ifail;
float  alpha_single;
float  beta_single;
double  alpha;
double  beta;
double  wkopt;
double  *work;
double  vl;
double  vu;

int idshead, *idsorder, *predorder, *usedids, *usedpreds, *indexer, *indexer2,
    *retain;
char idsfile[500] = "blank", **kinids, **kinids2, **kinids3, **wantids,
     **wantids2;
char  **wantpreds;
char  **wantpreds2;
char  **kinpreds;

int  readint;
int  readint2;
int  readint3;
int  *readints;
int  writeint;
float  *readfloats;
float  writefloat;
float  *writefloats;
double  readdouble;
double  readdouble2;
double  *readdoubles;
char  readchar;
char  readchar2;
char  *rs;
char **readstrings, readstring[500], readstring2[500], readstring3[500],
    readstring4[500], readstring5[500], readstring6[500];
char  **writestrings;
char  writestring[500];
char  cmd[2000];
char  cmd2[2000];

char filename[500], filename2[500], filename3[500], filename4[500],
    filename5[500], filename6[500], filename7[500], filename8[500],
    filename9[500], filename10[500];
FILE  *input;
FILE  *input2;
FILE  *input3;
FILE *output, *output2, *output3, *output4, *output5, *output6, *output7,
    *output8, *output9, *output10;

DIR  *dir;
gzFile  datainputgz;
struct  stat statstruct;

struct  sorting_double *dptrs;
struct  sorting_string *sptrs;

/////////////////////////// 

