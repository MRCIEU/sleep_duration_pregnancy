*input MoBa phenotype data
import delimited "\path\to\the\data", varnames(1) case(upper) stringcols(1 2 3) ///
numericcols(4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 ///
40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 ///
78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 ///
112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 ///
141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 ///
170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 ///
199 200 201 202 203) clear 

*clean sleep duration
*4, >=10 hour; 3, 8-9 hour; 2, 6-7 hour; 1, <=5 hour
gen sleep =.
replace sleep =4 if (CC1027==1)
replace sleep =3 if (CC1027==2)
replace sleep =2 if (CC1027==3)
replace sleep =1 if (CC1027==4)|(CC1027==5)|(CC1027==9)
tab sleep, mis

*clean covariates
replace MORS_ALDER =16 if (MORS_ALDER==917)
replace MORS_ALDER =46 if (MORS_ALDER==945)

*MI
keep if (sleep !=.)|(MORS_ALDER !=.)|(AA1124 !=.)|(AA1126 !=.)|(SMOKE !=.)|(SMOKEFTHR !=.)|(SMOKERPLC !=.)|(ALCOHOL !=.)| ///
(BMI !=.)|(INCOME !=.)|(VEKT !=.)|(PARITET_5 !=.)

mi set mlong
mi register imputed sleep MORS_ALDER AA1124 AA1126 SMOKE SMOKEFTHR ALCOHOL BMI INCOME VEKT PARITET_5
mi register regular SMOKERPLC
mi impute chained (mlogit) sleep (regress) MORS_ALDER (ologit) AA1124 (ologit) AA1126 (logit) SMOKE (logit) SMOKEFTHR (logit) ///
ALCOHOL (regress) BMI (ologit) INCOME (regress) VEKT (ologit) PARITET_5 = i.SMOKERPLC, dryrun
mi impute chained (mlogit) sleep (regress) MORS_ALDER (ologit) AA1124 (ologit) AA1126 (logit) SMOKE (logit) SMOKEFTHR (logit) ///
ALCOHOL (regress) BMI (ologit) INCOME (regress) VEKT (ologit) PARITET_5 = i.SMOKERPLC, add(100) rseed(1234)

*save imputed data
save "\path\to\the\data", replace

*analyze imputed data
use "\path\to\the\data", clear
*mi import mlong, automatic
mi describe

gen sleepL = .
replace sleepL=3.5 if (sleep==1)
replace sleepL=6.5 if (sleep==2)
replace sleepL=8.5 if (sleep==3)
replace sleepL=11 if (sleep==4)


mi estimate: regress VEKT sleepL MORS_ALDER i.AA1124 i.SMOKE i.ALCOHOL BMI i.INCOME PARITET_5
mi estimate: regress VEKT ib3.sleep MORS_ALDER i.AA1124 i.SMOKE i.ALCOHOL BMI i.INCOME PARITET_5


*likelihood ratio test
gen sleep3=sleep==3
gen sleep4=sleep==4
mi estimate: regress VEKT sleepL MORS_ALDER i.AA1124 i.SMOKE i.ALCOHOL BMI i.INCOME PARITET_5 sleep3 sleep4 
mi test sleep3 sleep4

