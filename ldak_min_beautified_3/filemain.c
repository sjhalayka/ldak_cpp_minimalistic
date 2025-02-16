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
// Functions for reading stuff
// general syntax - thing to read, objects to be filled, length, keeps, other
// stuff req to read, misc
///////////////////////////
void read_bimfile(char *bimfile, int *chr, char **preds, double *cm, double *bp,
                  char **along1, char **along2, char *al1, char *al2,
                  int length, int type, double window_cm, int flag, int exlong)
// type=0 - keep quiet, type=1 - warn about everything
// flag=0 - positions must be in order, flag=1 - do not
// exlong=0 - allowed long alleles, exlong=1 - not allowed
{
  int j, count, bcount, lcount, mcount, ncount, scount;
  char *rc, *rs, *rm, *rbp, *ra1, *ra2;
  FILE *input;
  rc = malloc(sizeof(char) * 10000000);
  rs = malloc(sizeof(char) * 10000000);
  rm = malloc(sizeof(char) * 10000000);
  rbp = malloc(sizeof(char) * 10000000);
  ra1 = malloc(sizeof(char) * 10000000);
  ra2 = malloc(sizeof(char) * 10000000);
  if ((input = fopen(bimfile, "r")) == NULL) {
    printf("Error opening %s\n\n", bimfile);
    exit(1);
  }
  bcount = 0;
  mcount = 0;
  ncount = 0;
  scount = 0;
  lcount = 0;
  for (j = 0;j < length;j++) {;
    if (fscanf(input, "%s %s %s %s %s %s ", rc, rs, rm, rbp, ra1, ra2) != 6) {
      printf("Error reading Row %d of %s\n\n", j + 1, bimfile);
      exit(1);
    }
    // save chr
    chr[j] = -1;
    if (strcmp(rc, "X") == 0) {
      chr[j] = 23;
    }
    if (strcmp(rc, "Y") == 0) {
      chr[j] = 24;
    }
    if (strcmp(rc, "XY") == 0) {
      chr[j] = 25;
    }
    if (strcmp(rc, "MT") == 0) {
      chr[j] = 26;
    }
    if (strcmp(rc, "0") == 0) {
      chr[j] = 0;
    }
    if (chr[j] == -1) // not found yet
    {
      chr[j] = atoi(rc);
      if (chr[j] <= 0) // so not valid
      {
        printf("Error, Predictor %s has chromosome %s  all values must be a "
               "positive integer, X (23), Y (24), XY (25), MT (26) or 0\n\n",
               rs, rc);
        exit(1);
      }
      if (chr[j] > 26 && type == 1) // unexpectedly large
      {
        if (bcount < 5) {
          printf("Warning, Predictor %s has chromosome %s (the largest "
                 "expected for humans is 26)\n",
                 rs, rc);
        }
        bcount++;
      }
    }
    if (j > 0 && flag == 0) // check chromosome not smaller than previous one
    {
      if (chr[j] < chr[j - 1]) {
        printf(
            "Error, chromosome for %s (%d) is lower than that for %s (%d)\n\n",
            rs, chr[j], preds[j - 1], chr[j - 1]);
        exit(1);
      }
    }
    // save name
    copy_string(preds, j, rs);
    // save cm and maybe warn if negative
    cm[j] = atof(rm);
    if (cm[j] < 0 && type == 1) {
      if (mcount < 5) {
        printf("Warning, Predictor %s has a negative genetic distance (%s)\n",
               rs, rm);
      }
      mcount++;
    }
    if (window_cm != -9999 && j > 0 &&
        flag == 0) // check cm consistent with previous one
    {
      if (chr[j] == chr[j - 1] && cm[j] < cm[j - 1]) {
        printf("Error, genetic distance for %s (%.2f) is lower than that for "
               "%s (%.2f)\n",
               rs, cm[j], preds[j - 1], cm[j - 1]);
        exit(1);
      }
    }
    // save bp and maybe warn if non-positive
    bp[j] = atof(rbp);
    if (bp[j] <= 0 && type == 1) {
      if (ncount < 5) {
        printf("Warning, Predictor %s has a non-positive basepair (%s) and "
               "will be ignored\n",
               rs, rbp);
      }
      ncount++;
    }
    if (j > 0 && bp[j] > 0 && flag == 0) // check bp consistent with previous
                                         // one and maybe warn if duplicate
    {
      if (chr[j] == chr[j - 1] && bp[j] < bp[j - 1]) {
        printf(
            "Error, basepair for %s (%.2f) is lower than that for %s (%.2f)\n",
            rs, bp[j], preds[j - 1], bp[j - 1]);
        exit(1);
      }
      if (chr[j] == chr[j - 1] && bp[j] == bp[j - 1] && type == 1) {
        if (scount < 5) {
          printf("Warning, %s and %s have the same basepair (%.2f)\n",
                 preds[j - 1], preds[j], bp[j]);
        }
        scount++;
      }
    }
    if (exlong == 0) // allowed long alleles - values of al1 and al2 redundant
    {
      // load long alleles
      copy_string(along1, j, ra1);
      copy_string(along2, j, ra2);
      // put first character into al1 and al2
      al1[j] = ra1[0];
      al2[j] = ra2[0];
    } else // not allowed long alleles - along1, along2, al1 and al2 redundant
           // if alleles long
    {
      if (strlen(ra1) + strlen(ra2) > 2) {
        if (type == 1) {
          if (lcount < 5) {
            printf("Warning, Predictor %s has multi-character alleles (%s and "
                   "%s) so will be ignored\n",
                   preds[j], ra1, ra2);
          }
          lcount++;
        }
        bp[j] = -1;
        strcpy(ra1, "X");
        strcpy(ra2, "Y");
      }
      copy_string(along1, j, ra1);
      copy_string(along2, j, ra2);
      al1[j] = ra1[0];
      al2[j] = ra2[0];
    }
  } // end of j loop
  fclose(input);
  if (window_cm !=
      -9999) // will be using maps, so check there are some valid values
  {
    count = 0;
    for (j = 0;j < length;j++) {;
      count += (cm[j] != 0);
    }
    if (count == 0) {
      printf("\nError, Column 3 of %s should contain genetic distances  either "
             "insert these or switch to physical distances (e.g,. replace "
             "\"--window-cm %.4f\" with \"--window-kb %.4f\")\n\n",
             bimfile, window_cm, window_cm * 1000);
      exit(1);
    }
  }
  if (bcount > 5) {
    printf("In total, %d predictors have chromosome values greater than 26\n",
           bcount);
  }
  if (mcount > 5) {
    printf("In total, %d predictors have negative genetic distances\n", mcount);
  }
  if (ncount > 5) {
    printf("In total, %d predictors have non-positive basepairs\n", ncount);
  }
  if (scount > 5) {
    printf("In total, %d pairs of predictors have the same basepair\n", scount);
  }
  if (lcount > 5) {
    printf("In total, %d predictors have multi-character alleles\n", lcount);
  }
  // if(bcount+mcount+ncount+scount+lcount>0){printf(" n");}
  free(rc);
  free(rs);
  free(rm);
  free(rbp);
  free(ra1);
  free(ra2);
} // end of read_bimfile
////////
void read_bimfile_bgen(char *bgenfile, int bgen_comp, int bgen_layout,
                       size_t *bgen_indexes, int *chr, char **preds, double *bp,
                       char **along1, char **along2, char *al1, char *al2,
                       int length, int type, int flag, int exlong)
// type=0 - keep quiet, type=1 - warn about everything
// flag=0 - positions must be in order, flag=1 - do not
// exlong=0 - allowed long alleles, exlong=1 - not allowed
{
  int j, k, count, count2, bcount, lcount, ncount, scount, tcount;
  char *rc, *rs, *ra1, *ra2, *rr;
  short idlen, readshort;
  int offset, readint, allen;
  unsigned int readuns;
  FILE *input;
  rc = malloc(sizeof(char) * 10000000);
  rs = malloc(sizeof(char) * 10000000);
  ra1 = malloc(sizeof(char) * 10000000);
  ra2 = malloc(sizeof(char) * 10000000);
  rr = malloc(sizeof(char) * 10000000);
  // open file
  if ((input = fopen(bgenfile, "rb")) == NULL) {
    printf("Error opening %s\n\n", bgenfile);
    exit(1);
  }
  // read offset, number of predictors and number of samples
  if (fread(&offset, 4, 1, input) != 1) {
    printf("Error reading first value of %s\n\n", bgenfile);
    exit(1);
  }
  fseeko(input, 8, SEEK_SET);
  if (fread(&count, 4, 1, input) != 1) {
    printf("Error reading third value of %s\n\n", bgenfile);
    exit(1);
  }
  if (fread(&count2, 4, 1, input) != 1) {
    printf("Error reading fourth value of %s\n\n", bgenfile);
    exit(1);
  }
  // skip to start of variants
  fseeko(input, 4 + offset, SEEK_SET);
  bcount = 0;
  ncount = 0;
  scount = 0;
  tcount = 0;
  lcount = 0;
  for (j = 0;j < length;j++) {;
    if ((j + 1) % 100000 == 0) {
      printf("Reading details for Predictor %d of %d\n", j + 1, length);
    }
    if (bgen_layout == 1) // read and check number of samples
    {
      if (fread(&readint, 4, 1, input) != 1) {
        printf("Error reading first value for Predictor %d of %s\n\n", j + 1,
               bgenfile);
        exit(1);
      }
      if (readint != count2) {
        printf("Error, %s seems to be corrupted (the numbers of samples are "
               "not consistent) - %d and %d\n\n",
               bgenfile, count2, readint);
        exit(1);
      }
    }
    // read short and string
    if (fread(&idlen, 2, 1, input) != 1) {
      printf("Error reading second value for Predictor %d of %s\n\n", j + 1,
             bgenfile);
      exit(1);
    }
    if (fread(rs, 1, idlen, input) != idlen) {
      printf("Error reading third value for Predictor %d of %s\n\n", j + 1,
             bgenfile);
      exit(1);
    }
    // read short and name
    if (fread(&idlen, 2, 1, input) != 1) {
      printf("Error reading fourth value for Predictor %d of %s\n\n", j + 1,
             bgenfile);
      exit(1);
    }
    if (fread(rs, 1, idlen, input) != idlen) {
      printf("Error reading fifth value for Predictor %d of %s\n\n", j + 1,
             bgenfile);
      exit(1);
    }
    rs[idlen] = '\0';
    // read short and chromosome
    if (fread(&idlen, 2, 1, input) != 1) {
      printf("Error reading sixth value for Predictor %d of %s\n\n", j + 1,
             bgenfile);
      exit(1);
    }
    if (fread(rc, 1, idlen, input) != idlen) {
      printf("Error reading seventh value for Predictor %d of %s\n\n", j + 1,
             bgenfile);
      exit(1);
    }
    rc[idlen] = '\0';
    // read position
    if (fread(&readuns, 4, 1, input) != 1) {
      printf("Error reading eighth value for Predictor %d of %s\n\n", j + 1,
             bgenfile);
      exit(1);
    }
    if (bgen_layout == 1) // have exactly two alleles
    {
      readshort = 2;
      // read first allele
      if (fread(&allen, 4, 1, input) != 1) {
        printf("Error reading first allele for Predictor %d of %s\n\n", j + 1,
               bgenfile);
        exit(1);
      }
      if (fread(ra1, 1, allen, input) != allen) {
        printf("Error reading first allele for Predictor %d of %s\n\n", j + 1,
               bgenfile);
        exit(1);
      }
      ra1[allen] = '\0';
      // read second allele
      if (fread(&allen, 4, 1, input) != 1) {
        printf("Error reading first allele for Predictor %d of %s\n\n", j + 1,
               bgenfile);
        exit(1);
      }
      if (fread(ra2, 1, allen, input) != allen) {
        printf("Error reading second allele for Predictor %d of %s\n\n", j + 1,
               bgenfile);
        exit(1);
      }
      ra2[allen] = '\0';
    } else // can be more (or less?) than two alleles
    {
      // read number of alleles
      if (fread(&readshort, 2, 1, input) != 1) {
        printf("Error reading ninth value for Predictor %d of %s\n\n", j + 1,
               bgenfile);
        exit(1);
      }
      if (readshort > 0) // read first allele
      {
        if (fread(&allen, 4, 1, input) != 1) {
          printf("Error reading first allele for Predictor %d of %s\n\n", j + 1,
                 bgenfile);
          exit(1);
        }
        if (fread(ra1, 1, allen, input) != allen) {
          printf("Error reading first allele for Predictor %d of %s\n\n", j + 1,
                 bgenfile);
          exit(1);
        }
        ra1[allen] = '\0';
      }
      if (readshort > 1) // read second allele
      {
        if (fread(&allen, 4, 1, input) != 1) {
          printf("Error reading second allele for Predictor %d of %s\n\n",
                 j + 1, bgenfile);
          exit(1);
        }
        if (fread(ra2, 1, allen, input) != allen) {
          printf("Error reading second allele for Predictor %d of %s\n\n",
                 j + 1, bgenfile);
          exit(1);
        }
        ra2[allen] = '\0';
      }
      for (k = 2;k < readshort;
           k++) // read the remaining alleles (will not use)
      {
        if (fread(&allen, 4, 1, input) != 1) {
          printf("Error reading second allele for Predictor %d of %s\n\n",
                 j + 1, bgenfile);
          exit(1);
        }
        if (fread(rr, 1, allen, input) != allen) {
          printf("Error reading Allele %d for Predictor %d of %s\n\n", k + 1,
                 j + 1, bgenfile);
          exit(1);
        }
      }
    }
    // save the file location
    bgen_indexes[j] = ftello(input);
    // skip to the next variant
    if (bgen_layout == 1) // for layout 1, length depends on compression
    {
      if (bgen_comp == 0) // must skip 6n
      {
        fseeko(input, 6 * count2, SEEK_CUR);
      } else // read length to skip
      {
        if (fread(&readint, 4, 1, input) != 1) {
          printf("Error reading storage size for Predictor %d of %s\n\n", j + 1,
                 bgenfile);
          exit(1);
        }
        fseeko(input, readint, SEEK_CUR);
      }
    } else // read length to skip (does not depend on compression)
    {
      if (fread(&readint, 4, 1, input) != 1) {
        printf("Error reading storage size for Predictor %d of %s\n\n", j + 1,
               bgenfile);
        exit(1);
      }
      fseeko(input, readint, SEEK_CUR);
    }
    // save chr
    chr[j] = -1;
    if (strcmp(rc, "X") == 0) {
      chr[j] = 23;
    }
    if (strcmp(rc, "Y") == 0) {
      chr[j] = 24;
    }
    if (strcmp(rc, "XY") == 0) {
      chr[j] = 25;
    }
    if (strcmp(rc, "MT") == 0) {
      chr[j] = 26;
    }
    if (strcmp(rc, "0") == 0) {
      chr[j] = 0;
    }
    if (chr[j] == -1) // not found yet
    {
      chr[j] = atoi(rc);
      if (chr[j] <= 0) // so not valid
      {
        printf("Error, Predictor %s has chromosome %s  all values must be a "
               "positive integer, X (23), Y (24), XY (25), MT (26) or 0\n\n",
               rs, rc);
        exit(1);
      }
      if (chr[j] > 26 && type == 1) // unexpectedly large
      {
        if (bcount < 5) {
          printf("Warning, Predictor %s has chromosome %s (the largest "
                 "expected for humans is 26)\n",
                 rs, rc);
        }
        bcount++;
      }
    }
    if (j > 0 && flag == 0) // check chromosome not smaller than previous one
    {
      if (chr[j] < chr[j - 1]) {
        printf(
            "Error, chromosome for %s (%d) is lower than that for %s (%d)\n\n",
            rs, chr[j], preds[j - 1], chr[j - 1]);
        exit(1);
      }
    }
    // save name
    copy_string(preds, j, rs);
    // save bp and maybe warn if non-positive
    bp[j] = (double)(readuns);
    if (bp[j] <= 0 && type == 1) {
      if (ncount < 5) {
        printf("Warning, Predictor %s has a non-positive basepair (%u) and "
               "will be ignored\n",
               rs, readuns);
      }
      ncount++;
    }
    if (j > 0 && bp[j] > 0 && flag == 0) // check bp consistent with previous
                                         // one and maybe warn if duplicate
    {
      if (chr[j] == chr[j - 1] && bp[j] < bp[j - 1]) {
        printf(
            "Error, basepair for %s (%.2f) is lower than that for %s (%.2f)\n",
            rs, bp[j], preds[j - 1], bp[j - 1]);
        exit(1);
      }
      if (chr[j] == chr[j - 1] && bp[j] == bp[j - 1] && type == 1) {
        if (scount < 5) {
          printf("Warning, %s and %s have the same basepair (%.2f)\n",
                 preds[j - 1], preds[j], bp[j]);
        }
        scount++;
      }
    }
    if (exlong == 0) // allowed long alleles - values of al1 and al2 redundant
    {
      // load long alleles
      copy_string(along1, j, ra1);
      copy_string(along2, j, ra2);
      // put first character into al1 and al2
      al1[j] = ra1[0];
      al2[j] = ra2[0];
    } else // not allowed long alleles - along1, along2, al1 and al2 redundant
           // if alleles long
    {
      if (strlen(ra1) + strlen(ra2) > 2) {
        if (type == 1) {
          if (lcount < 5) {
            printf("Warning, Predictor %s has multi-character alleles (%s and "
                   "%s) so will be ignored\n",
                   preds[j], ra1, ra2);
          }
          lcount++;
        }
        bp[j] = -1;
        strcpy(ra1, "X");
        strcpy(ra2, "Y");
      }
      copy_string(along1, j, ra1);
      copy_string(along2, j, ra2);
      al1[j] = ra1[0];
      al2[j] = ra2[0];
    }
  } // end of j loop
  // save total file size in final index
  bgen_indexes[length] = ftello(input);
  fclose(input);
  if (bcount > 5) {
    printf("In total, %d predictors have chromosome values greater than 26\n",
           bcount);
  }
  if (ncount > 5) {
    printf("In total, %d predictors have non-positive basepairs\n", ncount);
  }
  if (scount > 5) {
    printf("In total, %d pairs of predictors have the same basepair\n", scount);
  }
  if (tcount > 5) {
    printf("In total, %d predictors do not have exactly two alleles\n", lcount);
  }
  if (lcount > 5) {
    printf("In total, %d predictors have multi-character alleles\n", lcount);
  }
  // if(bcount+mcount+ncount+scount+lcount>0){printf(" n");}
  free(rc);
  free(rs);
  free(ra1);
  free(ra2);
  free(rr);
} // end of read_bimfile_bgen
///////////////////////////
int read_genfile(char *genfile, char **preds, double *bp, char **along1,
                 char **along2, char *al1, char *al2, int num_samples,
                 int genskip, int genheaders, int genprobs, int maxpreds,
                 int type, int flag, int exlong)
// type=0 - keep quiet, type=1 - warn about everything
// flag=0 - positions must be in order, flag=1 - do not
// exlong=0 - allowed long alleles, exlong=1 - not allowed
{
  int j, ncount, scount, lcount, size;
  char *rs, *rs2, *rbp, *ra1, *ra2;
  char *gzbuffer;
  gzFile inputgz;
  rs = malloc(sizeof(char) * 10000000);
  rs2 = malloc(sizeof(char) * 10000000);
  rbp = malloc(sizeof(char) * 10000000);
  ra1 = malloc(sizeof(char) * 10000000);
  ra2 = malloc(sizeof(char) * 10000000);
  if (genheaders < 4) {
    printf("Error B772D, genheader %d, please tell Doug\n\n", genheaders);
    exit(1);
  }
  // open, which checks size and skips header rows
  open_datagz(&inputgz, genfile, num_samples, genskip, genheaders, genprobs);
  // now can read
  size = 30000000 + num_samples * (genprobs * 20 + (genprobs == 0) * 4);
  gzbuffer = malloc(sizeof(char) * size);
  j = 0;
  ncount = 0;
  scount = 0;
  lcount = 0;
  while (gzgets(inputgz, gzbuffer, size) != NULL) {
    if (strlen(gzbuffer) == size - 1) {
      printf("Error reading %s  Row %d is longer (%d) than expected/allowed "
             "(%d)\nPlease tell Doug\n\n",
             genfile, genskip + 1, (int)strlen(gzbuffer), size);
      exit(1);
    }
    if ((genskip + j + 1) % 20000 == 0) {
      printf("Reading Row %d\n", genskip + j + 1);
    }
    if (genheaders == 4) {
      if (sscanf(gzbuffer, "%s %s %s %s ", rs, rbp, ra1, ra2) != 4) {
        printf("Error reading details on Row %d\n\n", genskip + j + 1);
        exit(1);
      }
    }
    if (genheaders == 5) {
      if (sscanf(gzbuffer, "%s %s %s %s %s ", rs2, rs, rbp, ra1, ra2) != 5) {
        printf("Error reading details on Row %d\n\n", genskip + j + 1);
        exit(1);
      }
    }
    if (genheaders == 6) {
      if (sscanf(gzbuffer, "%s %s %s %s %s %s ", rs2, rs, rs, rbp, ra1, ra2) !=
          6) {
        printf("Error reading details on Row %d\n\n", genskip + j + 1);
        exit(1);
      }
    }
    if (j < maxpreds) // saving
    {
      // save name
      copy_string(preds, j, rs);
      // save bp and maybe warn if non-positive
      bp[j] = atof(rbp);
      if (bp[j] <= 0 && type == 1) {
        if (ncount < 5) {
          printf("Warning, Predictor %s has a non-positive basepair (%s) and "
                 "will be ignored\n",
                 rs, rbp);
        }
        ncount++;
      }
      if (j > 0 && bp[j] > 0 && flag == 0) // check bp consistent with previous
                                           // one and maybe warn if duplicate
      {
        if (bp[j] < bp[j - 1]) {
          printf("Error, basepair for %s (%.2f) is lower than that for %s "
                 "(%.2f)\n",
                 rs, bp[j], preds[j - 1], bp[j - 1]);
          exit(1);
        }
        if (bp[j] == bp[j - 1] && type == 1) {
          if (scount < 5) {
            printf("Warning, %s and %s have the same basepair (%.2f)\n",
                   preds[j - 1], preds[j], bp[j]);
          }
          scount++;
        }
      }
      if (exlong == 0) // allowed long alleles - values of al1 and al2 redundant
      {
        // load long alleles
        copy_string(along1, j, ra1);
        copy_string(along2, j, ra2);
        // put first character into al1 and al2
        al1[j] = ra1[0];
        al2[j] = ra2[0];
      } else // not allowed long alleles - along1, along2, al1 and al2 redundant
             // if alleles long
      {
        if (strlen(ra1) + strlen(ra2) > 2) {
          if (type == 1) {
            if (lcount < 5) {
              printf("Warning, Predictor %s has multi-character alleles (%s "
                     "and %s) so will be ignored\n",
                     preds[j], ra1, ra2);
            }
            lcount++;
          }
          bp[j] = -1;
          strcpy(ra1, "X");
          strcpy(ra2, "Y");
        }
        copy_string(along1, j, ra1);
        copy_string(along2, j, ra2);
        al1[j] = ra1[0];
        al2[j] = ra2[0];
      }
    } // end of saving
    j++;
  } // end of while loop
  gzclose(inputgz);
  if (j <= maxpreds) {
    if (ncount > 5) {
      printf("In total, %d predictors have non-positive basepairs\n", ncount);
    }
    if (scount > 5) {
      printf("In total, %d pairs of predictors have the same basepair\n",
             scount);
    }
    if (lcount > 5) {
      printf("In total, %d predictors have multi-character alleles\n", lcount);
    }
  }
  free(rs);
  free(rs2);
  free(rbp);
  free(ra1);
  free(ra2);
  free(gzbuffer);
  return (j);
} // end of read_genfile
///////////////////////////
void read_bgen_headers(char *bgenfile, int *bgen_preds, int *bgen_samples,
                       int *bgen_comp, int *bgen_layout, int *bgen_ids) {
  char magic[4];
  unsigned char flags[4];
  int offset, hblock;
  FILE *input;
  if ((input = fopen(bgenfile, "rb")) == NULL) {
    printf("Error opening %s\n\n", bgenfile);
    exit(1);
  }
  fseeko(input, 0, SEEK_SET);
  // read offset and header length
  if (fread(&offset, 4, 1, input) != 1) {
    printf("Error reading first value of %s\n\n", bgenfile);
    exit(1);
  }
  if (fread(&hblock, 4, 1, input) != 1) {
    printf("Error reading second value of %s\n\n", bgenfile);
    exit(1);
  }
  // read number of preds and number of samples
  if (fread(bgen_preds, 4, 1, input) != 1) {
    printf("Error reading third value of %s\n\n", bgenfile);
    exit(1);
  }
  if (fread(bgen_samples, 4, 1, input) != 1) {
    printf("Error reading fourth value of %s\n\n", bgenfile);
    exit(1);
  }
  // check magic number, then skip free space
  if (fread(magic, 1, 4, input) != 4) {
    printf("Error reading fifth value of %s\n\n", bgenfile);
    exit(1);
  }
  if (magic[0] != 'b' && magic[0] != '0') {
    printf("Error reading %s  the first magic number should be b or 0 (not "
           "%c)\n\n",
           bgenfile, magic[0]);
    exit(1);
  }
  if (magic[1] != 'g' && magic[1] != '0') {
    printf("Error reading %s  the second magic number should be g or 0 (not "
           "%c)\n\n",
           bgenfile, magic[1]);
    exit(1);
  }
  if (magic[2] != 'e' && magic[2] != '0') {
    printf("Error reading %s  the third magic number should be e or 0 (not "
           "%c)\n\n",
           bgenfile, magic[2]);
    exit(1);
  }
  if (magic[3] != 'n' && magic[3] != '0') {
    printf("Error reading %s  the fourth magic number should be n or 0 (not "
           "%c)\n\n",
           bgenfile, magic[3]);
    exit(1);
  }
  fseeko(input, 4 + hblock - 4, SEEK_SET);
  // read the flags
  if (fread(flags, 1, 4, input) != 4) {
    printf("Error reading final header value from %s\n\n", bgenfile);
    exit(1);
  }
  *bgen_comp = flags[0] & 3;
  *bgen_layout = (flags[0] >> 2) & 3;
  *bgen_ids = (int)((flags[3] >> 7) & 1);
  if (*bgen_comp == 2) {
    printf("Sorry, %s uses ZSTD compression, but LDAK only allows for zlib "
           "compression  if you made the file using PLINK, please remake using "
           "\"--export bgen-1.2\" instead of \"--export bgen-1.3\" (if this is "
           "a problem, please ask Doug to add the feature\n\n",
           bgenfile);
    exit(1);
  }
  fclose(input);
} // end of read_bgen_headers
///////////////////////////
void check_respfile(char *respfile, int *usedids, int num_samples,
                    char **allids3, int num_resps_use, int *keepresps,
                    int num_resps, double missingvalue, int type)
// type=0 - no missing allowed, type=1 - missing allowed, but exclude samples
// missing all phenotypes
{
  int i, m, count, count2, count3, head;
  int *founds, *indexer, *indexer2;
  char **wantids;
  double *resptemp;
  char readchar, readstring[500];
  FILE *input;
  // want to tally how many responses present for each sample
  founds = malloc(sizeof(int) * num_samples);
  for (i = 0;i < num_samples;i++) {;
    founds[i] = 0;
  }
  // see which individuals are available
  head = check_head_ids(respfile, 0);
  count = countrows_plus(respfile, 2 + num_resps) - head;
  printf("Checking responses for %d samples from %s\n", count, respfile);
  wantids = malloc(sizeof(char *) * count);
  read_ids(respfile, NULL, NULL, wantids, count, NULL, head, 0);
  indexer = malloc(sizeof(int) * count);
  indexer2 = malloc(sizeof(int) * count);
  count2 = find_strings(wantids, count, allids3, num_samples, indexer, indexer2,
                        respfile, NULL, NULL, NULL, 3);
  if (count2 == 0) {
    printf("Error, can not find any of these samples\n\n");
    exit(1);
  }
  // open respfile and skip the header row if present
  if ((input = fopen(respfile, "r")) == NULL) {
    printf("Error opening %s\n\n", respfile);
    exit(1);
  }
  if (head == 1) {
    readchar = 0;
    while (readchar != 10) {
      readchar = 10;
      (void)fscanf(input, "%c", &readchar);
    }
  }
  resptemp = malloc(sizeof(double) * num_resps);
  count3 = 0;
  for (i = 0;i < count;i++) {;
    if (fscanf(input, "%s %s ", readstring, readstring) != 2) {
      printf("Error reading IDs on row %d of %s\n\n", i + 1, respfile);
      exit(1);
    }
    if (i == indexer[count3]) // using this row
    {
      for (m = 0;m < num_resps;m++) {;
        if (fscanf(input, "%s ", readstring) != 1) {
          printf("Error reading Element %d of Row %d of %s\n\n", 3 + m,
                 head + i + 1, respfile);
          exit(1);
        }
        if (strcmp(readstring, "NA") == 0) {
          resptemp[m] = missingvalue;
        } else {
          if (sscanf(readstring, "%lf%c", resptemp + m, &readchar) != 1) {
            printf("Error reading %s  Element %d of Row %d does not appear to "
                   "be numeric (%s)\n\n",
                   respfile, 3 + m, head + i + 1, readstring);
            exit(1);
          }
        }
      }
      for (m = 0;m < num_resps_use;m++) {;
        founds[indexer2[count3]] += (resptemp[keepresps[m]] != missingvalue);
      }
      count3++;
      if (count3 == count2) {
        break;
      }
    } else // not using this row
    {
      readchar = 0;
      while (readchar != 10) {
        readchar = 10;
        (void)fscanf(input, "%c", &readchar);
      }
    }
  }
  fclose(input);
  for (i = 0;i < count;i++) {;
    free(wantids[i]);
  }
  free(wantids);
  free(indexer);
  free(indexer2);
  free(resptemp);
  if (type == 0) // see if valid
  {
    for (i = 0;i < num_samples;i++) {;
      if (usedids[i] == 1 && founds[i] != 0 &&
          founds[i] != num_resps_use) // so must have num_resps_use>1
      {
        printf("Error, to analyse multiple phenotypes, each sample must have "
               "either all phenotypes present or all missing  you should "
               "instead either analyse each phenotype separately, or use "
               "\"--dentist YES\" to pad missing values\n\n");
        exit(1);
      }
    }
  }
  // remove samples with no phenotypes
  count = 0;
  for (i = 0;i < num_samples;i++) {;
    count += usedids[i];
  }
  count2 = 0;
  for (i = 0;i < num_samples;i++) {;
    usedids[i] = (usedids[i] == 1 && founds[i] > 0);
    count2 += usedids[i];
  }
  if (count2 == 0) {
    printf("Error, all phenotypes are missing\n\n");
    exit(1);
  }
  if (count2 < count) {
    printf("Due to missing phenotypic values, the number of samples is reduced "
           "from %d to %d\n",
           count, count2);
  }
  printf("\n");
  free(founds);
} // end of check_respfile
////////
void read_respfile(char *respfile, double *resp, int ns, char **ids3,
                   int num_resps_use, int *keepresps, int num_resps,
                   int *respcounts, double missingvalue, double binary,
                   int type)
// type=0 - reading only one phenotye, type=1/2 - reading first/second of a pair
// of phenotypes
{
  int i, m, count, count2, count3, indcount, head;
  int n0, n1, n2, nm, nr;
  double sum, sumsq, mean, var;
  int *indexer, *indexer2;
  char **wantids;
  double *resptemp;
  char readchar, readstring[500];
  FILE *input;
  // see which individuals are available
  head = check_head_ids(respfile, 0);
  count = countrows_plus(respfile, 2 + num_resps) - head;
  if (type == 0) {
    printf("Reading phenotypes for %d samples from %s\n", count, respfile);
  }
  if (type == 1) {
    printf("Reading first phenotype for %d samples from %s\n", count, respfile);
  }
  if (type == 2) {
    printf("Reading second phenotype for %d samples from %s\n", count,
           respfile);
  }
  wantids = malloc(sizeof(char *) * count);
  read_ids(respfile, NULL, NULL, wantids, count, NULL, head, 0);
  indexer = malloc(sizeof(int) * count);
  indexer2 = malloc(sizeof(int) * count);
  count2 = find_strings(wantids, count, ids3, ns, indexer, indexer2, respfile,
                        NULL, NULL, NULL, 3);
  if (count2 == 0) {
    printf("Error, can not find any of these samples\n\n");
    exit(1);
  }
  // set responses to missing
  for (m = 0;m < num_resps_use;m++) {;
    for (i = 0;i < ns;i++) {;
      resp[i + m * ns] = missingvalue;
    }
  }
  // open respfile and skip the header row if present
  if ((input = fopen(respfile, "r")) == NULL) {
    printf("Error opening %s\n\n", respfile);
    exit(1);
  }
  if (head == 1) {
    readchar = 0;
    while (readchar != 10) {
      readchar = 10;
      (void)fscanf(input, "%c", &readchar);
    }
  }
  resptemp = malloc(sizeof(double) * num_resps);
  count3 = 0;
  for (i = 0;i < count;i++) {;
    if (fscanf(input, "%s %s ", readstring, readstring) != 2) {
      printf("Error reading IDs on row %d of %s\n\n", i + 1, respfile);
      exit(1);
    }
    if (i == indexer[count3]) // using this row
    {
      for (m = 0;m < num_resps;m++) {;
        if (fscanf(input, "%s ", readstring) != 1) {
          printf("Error reading Element %d of Row %d of %s\n\n", 3 + m,
                 head + i + 1, respfile);
          exit(1);
        }
        if (strcmp(readstring, "NA") == 0) {
          resptemp[m] = missingvalue;
        } else {
          if (sscanf(readstring, "%lf%c", resptemp + m, &readchar) != 1) {
            printf("Error reading %s  Element %d of Row %d does not appear to "
                   "be numeric (%s)\n\n",
                   respfile, 3 + m, head + i + 1, readstring);
            exit(1);
          }
        }
      }
      for (m = 0;m < num_resps_use;m++) {;
        resp[indexer2[count3] + m * ns] = resptemp[keepresps[m]];
      }
      count3++;
      if (count3 == count2) {
        break;
      }
    } else // not using this row
    {
      readchar = 0;
      while (readchar != 10) {
        readchar = 10;
        (void)fscanf(input, "%c", &readchar);
      }
    }
  }
  fclose(input);
  for (i = 0;i < count;i++) {;
    free(wantids[i]);
  }
  free(wantids);
  free(indexer);
  free(indexer2);
  free(resptemp);
  ////////
  // check whether binary
  for (m = 0;m < num_resps_use;m++) {;
    n0 = 0;
    n1 = 0;
    n2 = 0;
    nm = 0;
    for (i = 0;i < ns;i++) {;
      if (resp[i + m * ns] == 0) {
        n0++;
      }
      if (resp[i + m * ns] == 1) {
        n1++;
      }
      if (resp[i + m * ns] == 2) {
        n2++;
      }
      if (resp[i + m * ns] == missingvalue) {
        nm++;
      }
    }
    nr = ns - n0 - n1 - n2 - nm;
    if (binary == 1 && nr > 0) {
      printf("Error, Phenotype %d is not binary\n\n", keepresps[m] + 1);
      exit(1);
    }
    if (binary == 1 && n0 > 0 && n1 > 0 && n2 > 0) {
      printf("Error, Phenotype %d takes three different values (note that LDAK "
             "does NOT treat 0 as missing)\n\n",
             keepresps[m] + 1);
      exit(1);
    }
    if (nr == 0 && n0 > 0 && n1 > 0 && n2 == 0) // cases=1, controls=0
    {
      if (m < 10) {
        printf("Response %d in %s is binary, with %d cases and %d controls",
               keepresps[m] + 1, respfile, n1, n0);
        if (nm > 0) {
          printf(", %d missing", nm);
        }
        printf("\n");
      }
    }
    if (nr == 0 && n0 == 0 && n1 > 0 &&
        n2 > 0) // cases=2, controls=1 - subtract one to get cases=1, controls=0
    {
      if (m < 10) {
        printf("Response %d in %s is binary, with %d cases and %d controls",
               keepresps[m] + 1, respfile, n2, n1);
        if (nm > 0) {
          printf(" (and %d missing)", nm);
        }
        printf("\n");
      }
      for (i = 0;i < ns;i++) {;
        if (resp[i + m * ns] != missingvalue) {
          resp[i + m * ns]--;
        }
      }
    }
    if (nr == 0 && n0 > 0 && n1 == 0 &&
        n2 >
            0) // cases=2, controls=0 - divide by two to get cases=1, controls=0
    {
      if (m < 10) {
        printf("Response %d in %s is binary, with %d cases and %d controls",
               keepresps[m] + 1, respfile, n2, n0);
        if (nm > 0) {
          printf(" (and %d missing)", nm);
        }
        printf("\n");
      }
      for (i = 0;i < ns;i++) {;
        if (resp[i + m * ns] != missingvalue) {
          resp[i + m * ns] *= 0.5;
        }
      }
    }
  } // end of m loop
  // get counts and check variances
  count = 0;
  for (m = 0;m < num_resps_use;m++) {;
    sum = 0;
    sumsq = 0;
    indcount = 0;
    for (i = 0;i < ns;i++) {;
      if (resp[i + m * ns] != missingvalue) {
        sum += resp[i + m * ns];
        sumsq += pow(resp[i + m * ns], 2);
        indcount++;
      }
    }
    if (indcount < 3) {
      printf(
          "Error reading %s  Phenotype %d has only %d non-missing values\n\n",
          respfile, keepresps[m] + 1, indcount);
      exit(1);
    }
    respcounts[m] = indcount;
    mean = sum / indcount;
    var = sumsq / indcount - pow(mean, 2);
    if (var == 0) {
      printf("Error reading %s  Phenotype %d has variance 0\n\n", respfile,
             keepresps[m] + 1);
      exit(1);
    }
    if (indcount < ns) // some missing values
    {
      if (m < 10) {
        printf("There are %d missing values values for Phenotype %d \n",
               ns - indcount, m + 1);
      }
      count++;
    }
  } // end of m loop
  if (num_resps_use > 1 && count > 0) {
    printf("In total, %d of the %d phenotypes have missing values\n", count,
           num_resps_use);
  }
  printf("\n");
} // end of read_respfile
///////////////////////////

