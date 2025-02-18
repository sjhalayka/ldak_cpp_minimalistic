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

// Data operations 

/////////////////////////// 

#include <memory>
using namespace std;

void open_datagz(gzFile* inputgz, char* datafile, int num_samples, int genskip,
    int genheaders, int genprobs) {
    int  j;
    int  count;
    int  size;
    int  size2;
    char  readchar;
    char  readchar2;
    char* gzbuffer;

    size = 30000000 + num_samples * (genprobs * 20 + (genprobs == 0) * 4);
    gzbuffer = malloc(sizeof(char) * size);

    // check size - open and read in until newline (or gzbuffer read) 
    if ((*inputgz = gzopen(datafile, "r")) == NULL) {
        printf("Error opening %s\n\n", datafile);
        exit(1);
    }

    if (gzgets(*inputgz, gzbuffer, size) == NULL) {
        printf("Error reading %s: file appears to be empty\n\n", datafile);
        exit(1);
    }
    size2 = strlen(gzbuffer);

    if (gzbuffer[0] == 9 || gzbuffer[0] == 10 || gzbuffer[0] == 32) {
        printf("Error, %s starts with a space or empty row\n", datafile);
        exit(1);
    }
    if (size2 == size - 1) {
        printf("Error reading %s: Row 1 is longer (%d) than expected/allowed "
            "(%d)\nPlease tell Doug\n\n",
            datafile, size2, size);
        exit(1);
    }

    // get number of elements in row 
    count = 1;
    readchar2 = gzbuffer[0];
    for (j = 1; j < size2 - 1; j++) {
        ;
        readchar = gzbuffer[j];
        if ((readchar2 == 9 || readchar2 == 32) && readchar != 9 &&
            readchar != 32) {
            count++;
        }
        readchar2 = readchar;
    }
    gzclose(*inputgz);

    if (genprobs > 0 && count != genheaders + num_samples * genprobs) {
        printf("Error reading %s: should contain %d columns, not %d (you may be "
            "able to fix this using \"--gen-headers\")\n\n",
            datafile, genheaders + num_samples * genprobs, count);
        exit(1);
    }

    if (genprobs == 0 && count != genheaders + num_samples * 2) {
        printf("Error reading %s: should contain %d columns, not %d (you may be "
            "able to fix this using \"--gen-headers\")\n\n",
            datafile, genheaders + num_samples * 2, count);
        exit(1);
    }

    // open and skip 
    if ((*inputgz = gzopen(datafile, "r")) == NULL) {
        printf("Error opening %s\n\n", datafile);
        exit(1);
    }
    count = 0;
    while (count < genskip) {
        if (gzgets(*inputgz, gzbuffer, size) == NULL) {
            printf("Error reading Header Row %d of %s\n\n", count + 1, datafile);
            exit(1);
        }
        size2 = strlen(gzbuffer);

        if (size2 >= size - 1) {
            printf("Error reading %s: Row %d is longer (%d) than expected/allowed "
                "(%d)\nPlease tell Doug\n\n",
                datafile, count + 1, size2, size);
            exit(1);
        }
        printf("Skipping Row %d of %s, I hope this is a header\n(Row begins: %c",
            count + 1, datafile, gzbuffer[0]);
        for (j = 1; j < size2 - 1; j++) {
            ;
            if (j < 50) {
                printf("%c", gzbuffer[j]);
            }
        }
        printf(")\n");
        count++;
    }
    if (count > 0) {
        printf("\n");
    }

    free(gzbuffer);
} // end of open_datagz

/////////////////////////// 

void read_bed_wrapper(char* bedfile, double* data, double* centres,
    double* mults, double* sqdevs, double* rates,
    double* infos, int num_samples_use, int* keepsamps,
    int length, int* keeppreds, int num_samples,
    int num_preds, double missingvalue, int* bedzeros,
    int* bedones, int* bedtwos, int type, int maxthreads) {
    int  thread;
    int  threadstart;
    int  threadend;
    int  threadlength;

    threadlength = (length - 1) / maxthreads + 1;

#pragma omp parallel for private(thread, threadstart, threadend)  schedule(static, 1)
    for (thread = 0; thread < maxthreads; thread++) {
        ;
        threadstart = thread * threadlength;
        threadend = (thread + 1) * threadlength;
        if (threadend > length) {
            threadend = length;
        }

        read_bed_fast(bedfile, data + (size_t)threadstart * num_samples_use,
            centres + threadstart, mults + threadstart,
            sqdevs + threadstart, rates + threadstart,
            infos + threadstart, num_samples_use, keepsamps,
            threadend - threadstart, keeppreds + threadstart, num_samples,
            num_preds, missingvalue, bedzeros, bedones, bedtwos, type);
    }
}

/////////////////////////// 

int read_data_fly(char* datafile, int dtype, double* data, float** probs,
    int num_samples_use, int* keepsamps, int start, int end,
    int* keeppreds_use, gzFile inputgz, size_t current,
    int num_samples, int num_preds, int genskip, int genheaders,
    int genprobs, size_t* bgen_indexes, double missingvalue,
    double threshold, double minprob, int nonsnp,
    int maxthreads) 
{
    int  thread;
    int  threadstart;
    int  threadend;
    int  threadlength;
 //   float*** threadprobs;

    vector<float> threadprobs_;

    //threadprobs_.resize(maxthreads);

    //threadprobs_.resize(threadprobs_.size() * 2);
    // 
    //auto ptr2 = make_unique<float*>(threadprobs_.data());
    ////auto ptr3 = make_unique<unique_ptr<float*>>(move(ptr2));

    //// These will cause segfault if probs is NULL
    //ptr2.get()[0] = probs[0] + (size_t)(threadstart - start) * num_samples_use;
    //ptr2.get()[1] = probs[1] + (size_t)(threadstart - start) * num_samples_use;

    //read_bgen_fly(
    //    datafile, data + (size_t)(threadstart - start) * num_samples_use,
    //    ptr2.get(), num_samples_use, keepsamps, threadstart,
    //    threadend, keeppreds_use, num_samples, num_preds, bgen_indexes,
    //    missingvalue, threshold, minprob);

    //return 0;






    if (dtype == 1 || dtype == 2 || dtype == 3 ||
        dtype == 4) // can read in parallel
    {
        threadlength = (end - start - 1) / maxthreads + 1;
        if (dtype == 2) {

            //threadprobs = malloc(sizeof(float**) * maxthreads);
        
            threadprobs_.resize(maxthreads);
        }

#pragma omp parallel for private(thread, threadstart, threadend)   schedule(static, 1)


        for (thread = 0; thread < maxthreads; thread++) {
            ;
            threadstart = start + thread * threadlength;
            threadend = start + (thread + 1) * threadlength;
            if (threadend > end) {
                threadend = end;
            }

            if (dtype == 1) {
                read_bed_fly(
                    datafile, data + (size_t)(threadstart - start) * num_samples_use,
                    num_samples_use, keepsamps, threadend - threadstart,
                    keeppreds_use + threadstart, num_samples, num_preds, missingvalue);
            }

            if (dtype == 2) {
                if (probs == NULL) {
                    read_bgen_fly(datafile,
                        data + (size_t)(threadstart - start) * num_samples_use,
                        NULL, num_samples_use, keepsamps, threadstart,
                        threadend, keeppreds_use, num_samples, num_preds,
                        bgen_indexes, missingvalue, threshold, minprob);
                }
                else {
                    //threadprobs[thread] = malloc(sizeof(float*) * 2);
                   
                    //size_t index0 = thread * threadprobs_.size() + 0;
                    //size_t index1 = thread * threadprobs_.size() + 1;

                    threadprobs_.resize(threadprobs_.size() * 2);

                    //threadprobs_[thread][0] =
                    //    probs[0] + (size_t)(threadstart - start) * num_samples_use;
                    //
                    //threadprobs_[thread][1] =
                    //    probs[1] + (size_t)(threadstart - start) * num_samples_use;

                    //const float* ptr1 = &threadprobs_[0];// +column_number * size_x * size_y];
                    //const float** ptr2 = &ptr1;

                   

                    auto ptr2 = make_unique<float*>(threadprobs_.data());
                    //auto ptr3 = make_unique<unique_ptr<float*>>(move(ptr2));

                    ptr2.get()[0] = probs[0] + (size_t)(threadstart - start) * num_samples_use;
                    ptr2.get()[1] = probs[1] + (size_t)(threadstart - start) * num_samples_use;

                    //read_bgen_fly(
                    //    datafile, data + (size_t)(threadstart - start) * num_samples_use,
                    //    threadprobs_[index0], num_samples_use, keepsamps, threadstart,
                    //    threadend, keeppreds_use, num_samples, num_preds, bgen_indexes,
                    //    missingvalue, threshold, minprob);

                    read_bgen_fly(
                        datafile, data + (size_t)(threadstart - start) * num_samples_use,
                        ptr2.get(), num_samples_use, keepsamps, threadstart,
                        threadend, keeppreds_use, num_samples, num_preds, bgen_indexes,
                        missingvalue, threshold, minprob);
                    
                    //free(threadprobs[thread]);
                }
            }

            if (dtype == 3) {
                read_sped_fly(
                    datafile, data + (size_t)(threadstart - start) * num_samples_use,
                    num_samples_use, keepsamps, threadstart, threadend, keeppreds_use,
                    num_samples, num_preds, missingvalue, threshold, nonsnp);
            }

            if (dtype == 4) {
                read_speed_fly(
                    datafile, data + (size_t)(threadstart - start) * num_samples_use,
                    num_samples_use, keepsamps, threadstart, threadend, keeppreds_use,
                    num_samples, num_preds, missingvalue, threshold, nonsnp);
            }
        }

        if (dtype == 2) {
           // free(threadprobs);
        }
    }

    if (dtype == 5) {
        (void)read_gen_fly(datafile, data, probs, num_samples_use, keepsamps, start,
            end, keeppreds_use, inputgz, current, num_samples,
            num_preds, genskip, genheaders, genprobs, missingvalue,
            threshold, minprob, nonsnp);
    }

    return (current + end);
}

/////////////////////////// 

void stand_data(double* data, double* centres, double* mults, double* sqdevs,
    double* rates, double* infos, int ns, int length,
    double missingvalue, double power, int gotcent, int hwestand,
    double* weights, int type)
    // type=0 - just calc stats (centres, mults, sqdevs, rates and infos) 
    // type=1 - calc stats, and standardize 
    // type=2 - standardize using provided centre and mult 
    // type=3 - calc stats, then set missing to mean 
{
    int  i;
    int  j;
    int  indcount;
    double  sum;
    double  sumsq;
    double  mean;
    double  var;
    double  value;

#pragma omp parallel for private(j, sum, sumsq, indcount, i, mean, var, value) schedule(static)



    for (j = 0; j < length; j++) {
        ;
        if (type == 0 || type == 1 ||
            type == 3) // calculate centre, mult, sqdev and rate
        {
            sum = 0;
            sumsq = 0;
            indcount = 0;
            for (i = 0; i < ns; i++) {
                ;
                if (data[(size_t)j * ns + i] != missingvalue) {
                    sum += data[(size_t)j * ns + i];
                    sumsq += pow(data[(size_t)j * ns + i], 2);
                    indcount++;
                }
            }

            if (indcount > 0) {
                mean = sum / indcount;
                var = sumsq / indcount - pow(mean, 2);
            }
            else {
                mean = 0;
                var = 0;
            }

            if (gotcent == 0) {
                centres[j] = mean;
            }
            else {
                mean = centres[j];
            }
            if (var > 0) {
                if (power == 0) {
                    mults[j] = 1;
                }
                else {
                    if (hwestand == 1) {
                        mults[j] = pow(mean * (1 - mean / 2), power / 2);
                    }
                    else {
                        mults[j] = pow(var * indcount / ns, power / 2);
                    }
                }
            }
            else {
                mults[j] = -9999;
            }
            sqdevs[j] = var * indcount / ns;
            rates[j] = (double)indcount / ns;
            value = centres[j] * (1 - centres[j] / 2);
            if (value > 0) {
                infos[j] = sqdevs[j] / value;
            }
            else {
                infos[j] = 0;
            }

            if (type == 3 && indcount < ns) // set missing to mean
            {
                for (i = 0; i < ns; i++) {
                    ;
                    if (data[(size_t)j * ns + i] == missingvalue) {
                        data[(size_t)j * ns + i] = mean;
                    }
                }
            }
        } // end of calculating centre, mult, sqdevs and rates

        if (type == 1 || type == 2) // standardize (with weights if provided)
        {
            if (weights != NULL) // using weights
            {
                if (mults[j] != -9999 && weights[j] > 0) // not trivial
                {
                    value = mults[j] * pow(weights[j], .5);
                    for (i = 0; i < ns; i++) {
                        ;
                        if (data[(size_t)j * ns + i] != missingvalue) {
                            data[(size_t)j * ns + i] =
                                (data[(size_t)j * ns + i] - centres[j]) * value;
                        }
                        else {
                            data[(size_t)j * ns + i] = 0;
                        }
                    }
                }
                else // trivial
                {
                    for (i = 0; i < ns; i++) {
                        ;
                        data[(size_t)j * ns + i] = 0;
                    }
                }
            }
            else // not using weights
            {
                if (mults[j] != -9999) // not trivial
                {
                    for (i = 0; i < ns; i++) {
                        ;
                        if (data[(size_t)j * ns + i] != missingvalue) {
                            data[(size_t)j * ns + i] =
                                (data[(size_t)j * ns + i] - centres[j]) * mults[j];
                        }
                        else {
                            data[(size_t)j * ns + i] = 0;
                        }
                    }
                }
                else // trivial
                {
                    for (i = 0; i < ns; i++) {
                        ;
                        data[(size_t)j * ns + i] = 0;
                    }
                }
            }
        } // end of type=1/2
    }   // end of j loop
} // end of stand_data

/////////////////////////// 

int extraction(int* usedpreds, int length, char** preds, int* chr,
    int* predorder, char* bpredfile, char* cpredfile, int onechr,
    char* onesnp, char* filename) {
    int  j;
    int  count;
    int  count2;
    int* indexer;
    char** wantpreds;

    count = 0;
    if (strcmp(bpredfile, "blank") != 0) {
        count = countrows(bpredfile);
        printf("Reading list of %d predictors to extract from %s\n", count,
            bpredfile);
        wantpreds = malloc(sizeof(char*) * count);
        read_strings(bpredfile, wantpreds, count, NULL, 1, 0);
    }
    if (strcmp(onesnp, "blank") != 0) {
        printf("Will only use Predictor %s\n", onesnp);
        count = 1;
        wantpreds = malloc(sizeof(char*));
        wantpreds[0] = malloc(sizeof(char) * (strlen(onesnp) + 1));
        strcpy(wantpreds[0], onesnp);
    }

    if (count > 0) {
        indexer = malloc(sizeof(int) * length);
        count2 = find_strings(preds, length, wantpreds, count, indexer, NULL, NULL,
            NULL, predorder, NULL, 3);
        if (count2 == 0) {
            if (strcmp(bpredfile, "blank") != 0) {
                printf("Error, none of these are in %s\n\n", filename);
                exit(1);
            }
            else {
                printf("Error, Predictor %s is not in %s\n\n", onesnp, filename);
                exit(1);
            }
        }

        if (count2 < count) {
            printf("Warning, only %d of these are in %s\n", count2, filename);
        }
        for (j = 0; j < count2; j++) {
            ;
            usedpreds[indexer[j]]++;
        }
        for (j = 0; j < length; j++) {
            ;
            usedpreds[j] = (usedpreds[j] == 2);
        }
        for (j = 0; j < count; j++) {
            ;
            free(wantpreds[j]);
        }
        free(wantpreds);
        free(indexer);
    }

    if (strcmp(cpredfile, "blank") != 0) // exclude predictors
    {
        count = countrows(cpredfile);
        printf("Reading list of %d predictors to exclude from %s (takes priority "
            "over any other predictor filtering)\n",
            count, cpredfile);
        wantpreds = malloc(sizeof(char*) * count);
        read_strings(cpredfile, wantpreds, count, NULL, 1, 0);

        indexer = malloc(sizeof(int) * length);
        count2 = find_strings(preds, length, wantpreds, count, indexer, NULL, NULL,
            NULL, predorder, NULL, 3);
        if (count2 == 0) {
            printf("Warning, none of these are in %s\n", filename);
        }
        if (count2 > 0 && count2 < count) {
            printf("Warning, only %d of these are in %s\n", count2, filename);
        }
        for (j = 0; j < count2; j++) {
            ;
            usedpreds[indexer[j]] = 0;
        }
        for (j = 0; j < count; j++) {
            ;
            free(wantpreds[j]);
        }
        free(wantpreds);
        free(indexer);
    }

    if (onechr > 0) {
        count = 0;
        for (j = 0; j < length; j++) {
            ;
            if (chr[j] == onechr) {
                count++;
            }
            else {
                usedpreds[j] = 0;
            }
        }
        if (count == 0) {
            printf("Error, there are no Chromosome %d predictors in %s\n\n", onechr,
                filename);
            exit(1);
        }
        printf("There are %d Chromosome %d predictors\n", count, onechr);
    }
    if (onechr == -1) {
        count = 0;
        for (j = 0; j < length; j++) {
            ;
            if (chr[j] > 0 && chr[j] < 23) {
                count++;
            }
            else {
                usedpreds[j] = 0;
            }
        }
        if (count == 0) {
            printf("Error, there are no autosomal predictors in %s\n\n", filename);
            exit(1);
        }
        printf("There are %d autosomal predictors\n", count);
    }
    if (onechr == -3) {
        count = 0;
        for (j = 0; j < length; j++) {
            ;
            if (chr[j] > 0 && chr[j] % 2 == 1) {
                count++;
            }
            else {
                usedpreds[j] = 0;
            }
        }
        if (count == 0) {
            printf("Error, there are no odd-chromosome predictors in %s\n\n",
                filename);
            exit(1);
        }
        printf("There are %d odd-chromosome predictors\n", count);
    }
    if (onechr == -2) {
        count = 0;
        for (j = 0; j < length; j++) {
            ;
            if (chr[j] > 0 && chr[j] % 2 == 0) {
                count++;
            }
            else {
                usedpreds[j] = 0;
            }
        }
        if (count == 0) {
            printf("Error, there are no even-chromosome predictors in %s\n\n",
                filename);
            exit(1);
        }
        printf("There are %d even-chromosome predictors\n", count);
    }

    count = 0;
    for (j = 0; j < length; j++) {
        ;
        count += usedpreds[j];
    }
    if (count == 0) {
        printf("Error, no predictors remain\n\n");
        exit(1);
    }

    printf("\n");

    return (count);
} // end of extraction

/////////////////////////// 

void change_coding(double* data, char* al1, char* al2, int ns, int length,
    int encoding, double missingvalue) {
    int  i;
    int  j;
    int  indcount;
    double  sum;

    char  readchar;

    for (j = 0; j < length; j++) {
        ;
        if (encoding == 2) // dominant - switch 1s for 2s
        {
            for (i = 0; i < ns; i++) {
                ;
                if (data[(size_t)j * ns + i] == 1) {
                    data[(size_t)j * ns + i] = 2;
                }
            }
        }
        if (encoding == 3) // recessive - switch 1s for 0s
        {
            for (i = 0; i < ns; i++) {
                ;
                if (data[(size_t)j * ns + i] == 1) {
                    data[(size_t)j * ns + i] = 0;
                }
            }
        }
        if (encoding == 4) // heterogeneous - 0/1/2 go to 0/2/0
        {
            for (i = 0; i < ns; i++) {
                ;
                if (data[(size_t)j * ns + i] != missingvalue) {
                    data[(size_t)j * ns + i] = 2 * (data[(size_t)j * ns + i] == 1);
                }
            }
            al1[j] = 'T';
            al2[j] = 'M';
        }
        if (encoding == 5) // ensure minor allele is a1 allele
        {
            sum = 0;
            indcount = 0;
            for (i = 0; i < ns; i++) {
                ;
                if (data[(size_t)j * ns + i] != missingvalue) {
                    sum += data[(size_t)j * ns + i];
                    indcount++;
                }
            }
            if (sum / indcount > 1) // switch
            {
                for (i = 0; i < ns; i++) {
                    ;
                    if (data[(size_t)j * ns + i] != missingvalue) {
                        data[(size_t)j * ns + i] = 2 - data[(size_t)j * ns + i];
                    }
                }
                readchar = al1[j];
                al1[j] = al2[j];
                al2[j] = readchar;
            }
        }
        if (encoding == 6) // missing or not missing
        {
            for (i = 0; i < ns; i++) {
                ;
                data[(size_t)j * ns + i] = (data[(size_t)j * ns + i] == missingvalue);
            }
            al1[j] = 'Y';
            al2[j] = 'N';
        }
    }
} // end of change_coding

/////////////////////////// 

