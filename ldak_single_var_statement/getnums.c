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
// Get names of samples / predictors and apply filterings - will not be here if
// merging be aware that find_strings can find multiple copies of stringb in
// string a not a problem for IDs, because no duplicates, but can be a problem
// for predictors
///////////////////////////
if (use_data == 1 || use_data == 2 || use_data == 3 || use_data == 4 ||
    num_kins > 0 || mode == 121 || mode == 123 || mode == 124 || mode == 129 ||
    mode == 130 || mode == 229 || mode == 230 ||
    mode == 177) // will be using sample ids (and will set idsfile)
{
  // set idsfile and get num_samples
  if (num_kins > 0 && mode != 122 && mode != 131 &&
      mode != 138) // start with intersect of kinships
  {
    sprintf(idsfile, "%s.grm.id", kinstems[0]);
    num_samples = total;
    if (num_kins == 1) {
      printf("Reading IDs for %d samples from %s\n\n", num_samples, idsfile);
    } else {
      printf("There are %d samples common to the %d kinship matrices\n\n",
             num_samples, num_kins);
    }
    count = countrows(idsfile);
    wantids = malloc(sizeof(char *) * count);
    read_ids(idsfile, NULL, NULL, wantids, count, NULL, 0, 0);
    indexer = malloc(sizeof(int) * count);
    count2 = find_strings(wantids, count, kinids, num_samples, indexer, NULL,
                          NULL, NULL, NULL, NULL, 1);
    if (count2 < num_samples) {
      printf("Error 78C  please tell Doug %d %d\n\n", count2, num_samples);
      exit(1);
    }
    allids1 = malloc(sizeof(char *) * num_samples);
    allids2 = malloc(sizeof(char *) * num_samples);
    allids3 = malloc(sizeof(char *) * num_samples);
    read_ids(filename, allids1, allids2, allids3, num_samples, indexer, 0, 0);
    for (i = 0;i < count;i++) {;
      free(wantids[i]);
    }
    free(wantids);
    free(indexer);
  } else {
    if (use_data == 1 || use_data == 2 || use_data == 3 ||
        use_data == 4) // start with famfile (or maybe bgenfile)
    {
      if (strcmp(famfile, datafile) != 0) {
        strcpy(idsfile, famfile);
        idshead = famhead;
        num_samples = countrows(idsfile) - idshead;
        printf("Reading IDs for %d samples from %s\n\n", num_samples, idsfile);
        allids1 = malloc(sizeof(char *) * num_samples);
        allids2 = malloc(sizeof(char *) * num_samples);
        allids3 = malloc(sizeof(char *) * num_samples);
        read_ids(idsfile, allids1, allids2, allids3, num_samples, NULL, idshead,
                 0);
      } else // read from  bgenfile
      {
        strcpy(idsfile, datafile);
        num_samples = bgen_samples;
        printf("Reading IDs for %d samples from %s\n\n", num_samples, datafile);
        allids1 = malloc(sizeof(char *) * num_samples);
        allids2 = malloc(sizeof(char *) * num_samples);
        allids3 = malloc(sizeof(char *) * num_samples);
        read_ids_bgen(datafile, allids1, allids2, allids3, num_samples, NULL,
                      1);
      }
    } else // start with response - must be null reml/he/pcgc, or no data/kins
           // cut-folds, or quant
    {
      strcpy(idsfile, respfile);
      idshead = check_head_ids(idsfile, 0);
      num_samples = countrows(idsfile) - idshead;
      printf("Reading IDs for %d samples from %s\n\n", num_samples, idsfile);
      allids1 = malloc(sizeof(char *) * num_samples);
      allids2 = malloc(sizeof(char *) * num_samples);
      allids3 = malloc(sizeof(char *) * num_samples);
      read_ids(idsfile, allids1, allids2, allids3, num_samples, NULL, idshead,
               0);
    }
  }
  // check no duplicates and get idsorder
  idsorder = malloc(sizeof(int) * num_samples);
  check_dups(allids3, num_samples, idsfile, idsorder, 1);
  ////////
  // usedids indicates which samples will be used
  usedids = malloc(sizeof(int) * num_samples);
  // start by assuming all num_samples are used
  for (i = 0;i < num_samples;i++) {;
    usedids[i] = 1;
  }
  if (strcmp(bsampfile, "blank") !=
      0) // keep samples (it is ok if bsampfile contains duplicates)
  {
    count = countrows(bsampfile);
    printf("Reading list of %d samples to keep from %s\n", count, bsampfile);
    wantids = malloc(sizeof(char *) * count);
    read_ids(bsampfile, NULL, NULL, wantids, count, NULL, 0, 0);
    indexer = malloc(sizeof(int) * num_samples);
    count2 = find_strings(allids3, num_samples, wantids, count, indexer, NULL,
                          NULL, NULL, idsorder, NULL, 3);
    if (count2 == 0) {
      printf("Error, none of these are in %s\n\n", idsfile);
      exit(1);
    }
    if (count2 < count) {
      printf("Warning, only %d of these are in %s\n", count2, idsfile);
    }
    printf("\n");
    for (i = 0;i < count2;i++) {;
      usedids[indexer[i]]++;
    }
    for (i = 0;i < num_samples;i++) {;
      usedids[i] = (usedids[i] == 2);
    }
    for (i = 0;i < count;i++) {;
      free(wantids[i]);
    }
    free(wantids);
    free(indexer);
  }
  if (num_subs > 0) // reduce to samples in at least one subset (it is ok if the
                    // subsets are not distinct)
  {
    count = 0;
    for (s = 0;s < num_subs;s++) {;
      sprintf(filename, "%s%d", subpref, s + 1);
      count += countrows(filename);
    }
    wantids = malloc(sizeof(char *) * count);
    count = 0;
    for (s = 0;s < num_subs;s++) {;
      sprintf(filename, "%s%d", subpref, s + 1);
      count2 = countrows(filename);
      read_ids(filename, NULL, NULL, wantids + count, count2, NULL, 0, 0);
      count += count2;
    }
    indexer = malloc(sizeof(int) * num_samples);
    count2 = find_strings(allids3, num_samples, wantids, count, indexer, NULL,
                          NULL, NULL, idsorder, NULL, 3);
    if (count2 == 0) {
      printf("Error, none of the samples in the %d subsets are in %s\n\n",
             num_subs, idsfile);
      exit(1);
    }
    for (i = 0;i < count2;i++) {;
      usedids[indexer[i]]++;
    }
    for (i = 0;i < num_samples;i++) {;
      usedids[i] = (usedids[i] == 2);
    }
    for (i = 0;i < count;i++) {;
      free(wantids[i]);
    }
    free(wantids);
    free(indexer);
  }
  if (strcmp(csampfile, "blank") !=
      0) // remove samples (it is ok if csampfile contains duplicates)
  {
    count = countrows(csampfile);
    printf("Reading list of %d samples to remove from %s", count, csampfile);
    if (strcmp(bsampfile, "blank") != 0) {
      printf(" (note that \"--remove\" takes priority over \"--keep\")");
    }
    printf("\n");
    wantids = malloc(sizeof(char *) * count);
    read_ids(csampfile, NULL, NULL, wantids, count, NULL, 0, 0);
    indexer = malloc(sizeof(int) * num_samples);
    count2 = find_strings(allids3, num_samples, wantids, count, indexer, NULL,
                          NULL, NULL, idsorder, NULL, 3);
    if (count2 < count) {
      if (count2 == 0) {
        printf("Warning, none of these are in %s\n", idsfile);
      } else {
        printf("Warning, only %d of these are in %s\n", count2, idsfile);
      }
    }
    printf("\n");
    for (i = 0;i < count2;i++) {;
      usedids[indexer[i]] = 0;
    }
    for (i = 0;i < count;i++) {;
      free(wantids[i]);
    }
    free(wantids);
    free(indexer);
  }
  count = 0;
  for (i = 0;i < num_samples;i++) {;
    count += (usedids[i] == 1);
  }
  if (count == 0) {
    printf("Error, after filtering, no samples remain\n\n");
    exit(1);
  }
  if (pad == 0) // exclude ind missing phenotypes (not for filter, blups,
                // reml-predict, family or scores, nor some linear)
  {
    if ((mode == 106 || mode == 107 || mode == 108 || mode == 109 ||
         mode == 112 || mode == 114 || mode == 116 || mode == 117 ||
         mode == 118 || mode == 119 || mode == 120 || mode == 121 ||
         mode == 123 || mode == 124 || mode == 126 || mode == 127 ||
         mode == 128 || (mode == 131 && trios == 0 && duos == 0) ||
         mode == 132 || mode == 133 || mode == 137 || mode == 138 ||
         mode == 140 || mode == 141 || mode == 145 || mode == 151 ||
         mode == 152 || mode == 153 || mode == 154 || mode == 156 ||
         mode == 160 || mode == 161 || mode == 163 || mode == 164 ||
         mode == 166 || mode == 167 || mode == 168 || mode == 169 ||
         mode == 170 || mode == 171 || mode == 177 || mode == 181 ||
         mode == 182 || mode == 183 || mode == 184 || mode == 185 ||
         mode == 186 || mode == 187 || mode == 188 || mode == 189 ||
         mode == 190 || mode == 191 || mode == 192 || mode == 193 ||
         mode == 194) &&
        strcmp(respfile, "blank") != 0) {
      check_respfile(respfile, usedids, num_samples, allids3, num_resps_use,
                     keepresps, num_resps, missingvalue, 0);
    }
  }
  if (pad == 2) // exclude ind without any phenotypes, then set pad=1
  {
    check_respfile(respfile, usedids, num_samples, allids3, num_resps_use,
                   keepresps, num_resps, missingvalue, 1);
    pad = 1;
  }
  ////////
  // see which samples remain, then do some mainly node-specific checks, then
  // load up ids
  num_samples_use = 0;
  for (i = 0;i < num_samples;i++) {;
    num_samples_use += (usedids[i] == 1);
  }
  if (num_samples_use < 3) {
    if (mode != 172) {
      printf("Error, unable to continue with fewer than three samples  come "
             "on, you can do better  )\n\n");
      exit(1);
    }
  }
  ids1 = malloc(sizeof(char *) * num_samples_use);
  ids2 = malloc(sizeof(char *) * num_samples_use);
  ids3 = malloc(sizeof(char *) * num_samples_use);
  count = 0;
  for (i = 0;i < num_samples;i++) {;
    if (usedids[i] == 1) {
      copy_string(ids1, count, allids1[i]);
      copy_string(ids2, count, allids2[i]);
      copy_string(ids3, count, allids3[i]);
      count++;
    }
  }
  for (i = 0;i < num_samples;i++) {;
    free(allids1[i]);
    free(allids2[i]);
    free(allids3[i]);
  }
  free(allids1);
  free(allids2);
  free(allids3);
  free(idsorder);
  free(usedids);
  ////////
  if (use_data == 1 || use_data == 2 ||
      use_data == 4) // set keepsamps, the indexes corresponding to ids3 (will
                     // be no duplicates)
  {
    if (strcmp(famfile, datafile) != 0) {
      num_samples = countrows(famfile) - famhead;
      wantids = malloc(sizeof(char *) * num_samples);
      read_ids(famfile, NULL, NULL, wantids, num_samples, NULL, famhead, 0);
    } else {
      num_samples = bgen_samples;
      wantids = malloc(sizeof(char *) * num_samples);
      read_ids_bgen(datafile, NULL, NULL, wantids, num_samples, NULL, 0);
    }
    keepsamps = malloc(sizeof(int) * num_samples_use);
    count = find_strings(wantids, num_samples, ids3, num_samples_use, keepsamps,
                         NULL, NULL, NULL, NULL, NULL, 3);
    if (count < num_samples_use) {
      printf("Error 3DE  please tell Doug %d %d\n\n", count, num_samples_use);
      exit(1);
    }
    for (i = 0;i < num_samples;i++) {;
      free(wantids[i]);
    }
    free(wantids);
  }
  if (use_data ==
      3) // using data only through toppreds - check some samples in common
  {
    if (strcmp(famfile, datafile) != 0) {
      num_samples = countrows(famfile) - famhead;
      wantids = malloc(sizeof(char *) * num_samples);
      read_ids(famfile, NULL, NULL, wantids, num_samples, NULL, famhead, 0);
    } else {
      num_samples = bgen_samples;
      wantids = malloc(sizeof(char *) * num_samples);
      read_ids_bgen(datafile, NULL, NULL, wantids, num_samples, NULL, 0);
    }
    count = find_strings(wantids, num_samples, ids3, num_samples_use, NULL,
                         NULL, NULL, NULL, NULL, NULL, 3);
    if (count == 0) {
      printf("Error, none of the samples are in %s\n\n", famfile);
      exit(1);
    }
    if (count < num_samples_use) {
      printf("Warning, only %d of the samples are in %s\n\n", count, famfile);
    }
    for (i = 0;i < num_samples;i++) {;
      free(wantids[i]);
    }
    free(wantids);
  }
  if (num_kins > 0) // free allocations from above
  {
    for (i = 0;i < total;i++) {;
      free(kinids[i]);
    }
    free(kinids);
  }
} // end of using sample ids
///////////////////////////
if (use_data == 1 || use_data == 2 || use_data == 3 ||
    use_data == 6) // will be using predictor ids, so set num_preds,
                   // num_preds_use and get predictor details
{
  // read in all predictors
  if (strcmp(bimfile, datafile) != 0) // have a separate bimfile
  {
    if (countcols(bimfile) != 6) {
      printf("Error, %s should have six columns (not %d)  columns should be "
             "chr, name, genetic distance, bp, A1, A2\n\n",
             bimfile, countcols(bimfile));
      exit(1);
    }
    num_preds = countrows_plus(bimfile, 6);
    printf("Reading details for %d predictors from %s\n", num_preds, bimfile);
    allchr = malloc(sizeof(int) * num_preds);
    allcm = malloc(sizeof(double) * num_preds);
    allbp = malloc(sizeof(double) * num_preds);
    allpreds = malloc(sizeof(char *) * num_preds);
    allalong1 = malloc(sizeof(char *) * num_preds);
    allalong2 = malloc(sizeof(char *) * num_preds);
    allal1 = malloc(sizeof(char) * num_preds);
    allal2 = malloc(sizeof(char) * num_preds);
    // currently always require predictors are ordered by position
    read_bimfile(bimfile, allchr, allpreds, allcm, allbp, allalong1, allalong2,
                 allal1, allal2, num_preds, 1, window_cm, 0, exlong);
  } else // no bimfile, must be dtypes 2 or 5
  {
    if (dtype == 2) // read from bgenfile (so can not have genetic distances)
    {
      num_preds = bgen_preds;
      printf("Reading details for %d predictors from %s (bgen layout %d, "
             "compression %d)\n",
             num_preds, bimfile, bgen_layout, bgen_comp);
      bgen_indexes = malloc(sizeof(size_t) * (num_preds + 1));
      allchr = malloc(sizeof(int) * num_preds);
      allcm = malloc(sizeof(double) * num_preds);
      allbp = malloc(sizeof(double) * num_preds);
      allpreds = malloc(sizeof(char *) * num_preds);
      allalong1 = malloc(sizeof(char *) * num_preds);
      allalong2 = malloc(sizeof(char *) * num_preds);
      allal1 = malloc(sizeof(char) * num_preds);
      allal2 = malloc(sizeof(char) * num_preds);
      for (j = 0;j < num_preds;j++) {;
        allcm[j] = 0;
      }
      // currently always require predictors are ordered by position
      read_bimfile_bgen(datafile, bgen_comp, bgen_layout, bgen_indexes, allchr,
                        allpreds, allbp, allalong1, allalong2, allal1, allal2,
                        num_preds, 1, 0, exlong);
    } else // read from genfile (so can not have genetic distances)
    {
      // assume at most maxpreds SNPs
      printf("Reading predictor details from %s\n", datafile);
      allchr = malloc(sizeof(int) * maxpreds);
      allcm = malloc(sizeof(double) * maxpreds);
      allbp = malloc(sizeof(double) * maxpreds);
      allpreds = malloc(sizeof(char *) * maxpreds);
      allalong1 = malloc(sizeof(char *) * maxpreds);
      allalong2 = malloc(sizeof(char *) * maxpreds);
      allal1 = malloc(sizeof(char) * maxpreds);
      allal2 = malloc(sizeof(char) * maxpreds);
      for (j = 0;j < maxpreds;j++) {;
        allchr[j] = oxchr;
        allcm[j] = 0;
      }
      // currently always require predictors are ordered by position
      num_preds = read_genfile(datafile, allpreds, allbp, allalong1, allalong2,
                               allal1, allal2, num_samples, genskip, genheaders,
                               genprobs, maxpreds, 1, 0, exlong);
      if (num_preds >
          maxpreds) // then file larger than allocated, will have to read again
      {
        for (j = 0;j < maxpreds;j++) {;
          free(allpreds[j]);
          free(allalong1[j]);
          free(allalong2[j]);
        }
        free(allpreds);
        free(allalong1);
        free(allalong2);
        free(allchr);
        free(allcm);
        free(allbp);
        free(allal1);
        free(allal2);
        allchr = malloc(sizeof(int) * num_preds);
        allcm = malloc(sizeof(double) * num_preds);
        allbp = malloc(sizeof(double) * num_preds);
        allpreds = malloc(sizeof(char *) * num_preds);
        allalong1 = malloc(sizeof(char *) * num_preds);
        allalong2 = malloc(sizeof(char *) * num_preds);
        allal1 = malloc(sizeof(char) * num_preds);
        allal2 = malloc(sizeof(char) * num_preds);
        for (j = 0;j < num_preds;j++) {;
          allchr[j] = oxchr;
          allcm[j] = 0;
        }
        (void)read_genfile(datafile, allpreds, allbp, allalong1, allalong2,
                           allal1, allal2, num_samples, genskip, genheaders,
                           genprobs, num_preds, 1, 0, exlong);
      }
      if (num_preds >= 20000) {
        printf("%s contains genotypes for %d SNPs\n", datafile, num_preds);
      }
    }
  }
  printf("\n");
  if (exlong == 0) // see if any long alleles
  {
    count = 0;
    for (j = 0;j < num_preds;j++) {;
      count +=
          (allbp[j] > 0 && (strlen(allalong1[j]) + strlen(allalong2[j]) > 2));
    }
    if (count > 0) {
      printf("Warning, %d predictors have multi-character alleles (you can "
             "remove these using \"--exclude-long-alleles YES\")\n\n",
             count);
    }
  }
  // check there are some predictors with positive bp
  count = 0;
  for (j = 0;j < num_preds;j++) {;
    count += (allbp[j] > 0);
  }
  if (count == 0) {
    if (exlong == 0) {
      printf("Error, there are no predictors with non-positive basepairs\n\n");
    } else {
      printf("Error, there are no predictors with non-positive basepairs and "
             "single-character alleles\n\n");
    }
    exit(1);
  }
  ////////
  // get predorder
  predorder = malloc(sizeof(int) * num_preds);
  (void)check_dups(allpreds, num_preds, bimfile, predorder, 0);
  // get to first retained predictor (must be at least one)
  found = 0;
  while (allbp[predorder[found]] <= 0) {
    found++;
  }
  // check if names are unique (after allowing for long alleles and negative
  // basepairs)
  wcount = 0;
  for (j = found + 1;j < num_preds;j++) {;
    if (allbp[predorder[j]] > 0) {
      if (strcmp(allpreds[predorder[j]], allpreds[predorder[found]]) ==
          0) // duplicate name
      {
        if (wcount < 5) {
          if (exsame == 0) {
            printf("Warning, at least two predictors are called %s (you can "
                   "remove these using \"--exclude-same-names YES\")\n",
                   allpreds[predorder[j]]);
          } else {
            printf("Warning, at least two predictors are called %s  these will "
                   "be ignored\n",
                   allpreds[predorder[j]]);
            allbp[predorder[found]] = -1;
            allbp[predorder[j]] = -1;
          }
        }
        wcount++;
        // if three or more with same name, don't want to keep warning
        while (j < num_preds - 1) {
          j++;
          if (allbp[predorder[j]] > 0) {
            if (strcmp(allpreds[predorder[j]], allpreds[predorder[found]]) !=
                0) {
              break;
            }
            if (exsame == 1) {
              allbp[predorder[j]] = -1;
            }
          }
        }
      }
      found = j;
    }
  }
  if (wcount > 5) {
    printf("In total, %d predictor names appear more than once\n", wcount);
  }
  if (wcount > 0) {
    printf("\n");
  }
  // check there are still some predictors with positive bp
  count = 0;
  for (j = 0;j < num_preds;j++) {;
    count += (allbp[j] > 0);
  }
  if (count == 0) {
    printf("Error, there are no predictors with unique names\n\n");
    exit(1);
  }
  if (count < num_preds) // update predorder so that missing predictors are not
                         // considered when using find_strings
  {
    count2 = 0;
    for (j = 0;j < num_preds;j++) {;
      if (allbp[predorder[j]] > 0) {
        if (count2 != j) {
          predorder[count2] = predorder[j];
        }
        count2++;
      }
    }
    for (j = count;j < num_preds;j++) {;
      predorder[j] = -1;
    }
  }
  if (strcmp(topfile, "blank") !=
      0) // check tops in the data (and were not excluded), and set tkeeppreds
  {
    wantpreds = malloc(sizeof(char *) * num_tops);
    read_strings(topfile, wantpreds, num_tops, NULL, 1, 0);
    tkeeppreds = malloc(sizeof(int) * num_preds);
    indexer = malloc(sizeof(int) * num_preds);
    usedpreds = malloc(sizeof(int) * num_tops);
    count = find_strings(allpreds, num_preds, wantpreds, num_tops, tkeeppreds,
                         indexer, NULL, topfile, predorder, NULL, 3);
    if (count == 0) {
      printf("Error, none of the %d top predictors are in the data\n\n",
             num_tops);
      exit(1);
    }
    // check all found (exactly once)
    for (j = 0;j < num_tops;j++) {;
      usedpreds[j] = 0;
    }
    for (j = 0;j < count;j++) {;
      usedpreds[indexer[j]]++;
    }
    for (j = 0;j < num_tops;j++) {;
      if (usedpreds[j] == 0) {
        printf("Error, the top predictor %s is not in the data\n\n",
               wantpreds[j]);
        exit(1);
      }
      if (usedpreds[j] > 1) {
        printf("Error, the top predictor %s appears %d times in the data\n\n",
               wantpreds[j], usedpreds[j]);
        exit(1);
      }
    }
    for (j = 0;j < num_tops;j++) {;
      free(wantpreds[j]);
    }
    free(wantpreds);
    free(indexer);
    free(usedpreds);
  }
  ////////
  // now set usedpreds
  usedpreds = malloc(sizeof(int) * num_preds);
  // blank predictors that have already been excluded
  for (j = 0;j < num_preds;j++) {;
    usedpreds[j] = (allbp[j] > 0);
  }
  if (extract == 1) // do predictor filtering
  {
    (void)extraction(usedpreds, num_preds, allpreds, allchr, predorder,
                     bpredfile, cpredfile, onechr, onesnp, bimfile);
  }
  ////////
  // set keeppreds
  keeppreds = malloc(sizeof(int) * num_preds);
  num_preds_use = 0;
  for (j = 0;j < num_preds;j++) {;
    if (usedpreds[j] == 1) {
      keeppreds[num_preds_use] = j;
      num_preds_use++;
    }
  }
  if (num_preds_use == 0) {
    printf("Error, after filtering predictors, none remain\n\n");
    exit(1);
  }
  // work out num_chr
  num_chr = 1;
  for (j = 1;j < num_preds_use;j++) {;
    if (allchr[keeppreds[j]] != allchr[keeppreds[j - 1]]) {
      num_chr++;
    }
  }
  if (strcmp(prsfile, "blank") != 0) {
    if (num_chr > num_chr2) {
      printf("Error, will be analyzing %d chromosomes, but only %d were "
             "included when making the PRS\n\n",
             num_chr, num_chr2);
      exit(1);
    }
  }
  // do some checks and set some values
  if (window_length == -1) {
    window_length = num_preds_use;
  }
  if (num_chr == 1 && (kvikstep == 2 || gctastep == 2 || faststep == 2) &&
      bychr == 1) // only one chromosome, so modify outfile
  {
    strcpy(filename, outfile);
    sprintf(outfile, "%s.chr%d", filename, allchr[keeppreds[0]]);
    printf("Warning, all predictors are on Chromosome %d, so will append "
           "\".chr%d\" to the output file stem (e.g., the main results will be "
           "saved in %s.assoc)  you can prevent this by adding \"--by-chr "
           "NO\"\n\n",
           allchr[keeppreds[0]], allchr[keeppreds[0]], outfile);
  }
  free(usedpreds);
} // end of use_data=1/2/3
///////////////////////////

