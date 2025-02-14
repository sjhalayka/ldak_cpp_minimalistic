/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//set data_length and keeppreds

///////////////////////////

//set data_length to num_preds_use (number after filtering), then overwrite (or update) when necessary
data_length=num_preds_use;
keeppreds_use=malloc(sizeof(int)*num_preds);
for(j=0;j<num_preds_use;j++){keeppreds_use[j]=keeppreds[j];}

////////

//now can reduce annotations and load cmbp (will contain either 1000000xcm or bp)
chr=malloc(sizeof(int)*data_length);
cm=malloc(sizeof(double)*data_length);
bp=malloc(sizeof(double)*data_length);
cmbp=malloc(sizeof(double)*data_length);
preds=malloc(sizeof(char*)*data_length);
along1=malloc(sizeof(char*)*data_length);
along2=malloc(sizeof(char*)*data_length);
al1=malloc(sizeof(char)*data_length);
al2=malloc(sizeof(char)*data_length);
for(j=0;j<data_length;j++)
{
chr[j]=allchr[keeppreds_use[j]];
cm[j]=allcm[keeppreds_use[j]];
bp[j]=allbp[keeppreds_use[j]];
if(window_cm!=-9999){cmbp[j]=1000000*allcm[keeppreds_use[j]];}
else{cmbp[j]=allbp[keeppreds_use[j]];}
copy_string(preds,j,allpreds[keeppreds_use[j]]);
copy_string(along1,j,allalong1[keeppreds_use[j]]);
copy_string(along2,j,allalong2[keeppreds_use[j]]);
al1[j]=allal1[keeppreds_use[j]];
al2[j]=allal2[keeppreds_use[j]];
}

///////////////////////////

