##############################################################################################################################################

## in addition, add chromosomal gains at early and late stages
tumors.with.17q <- rownames(sample.information.80x)[sample.information.80x$`17q.gain`=="T"]
tumors.with.whole.17.gain <- rownames(sample.information.80x)[sample.information.80x$`17.gain`=="T" & sample.information.80x$`17q.gain`=="F"]
tumors.with.7q <- rownames(sample.information.80x)[sample.information.80x$`7q.gain`=="T"]
tumors.with.whole.7.gain <- rownames(sample.information.80x)[sample.information.80x$`7.gain`=="T" & sample.information.80x$`7q.gain`=="F"]
tumors.with.1p.loss <- rownames(sample.information.80x)[sample.information.80x$`1p.deletion`=="T"]
tumors.with.1q.gain <- rownames(sample.information.80x)[sample.information.80x$`1q.gain`=="T"]
tumors.with.whole.1.gain <- rownames(sample.information.80x)[sample.information.80x$`1.gain`=="T" & sample.information.80x$`1q.gain`=="F"]
tumors.with.2p <- rownames(sample.information.80x)[sample.information.80x$`2p.gain`=="T"]
tumors.with.whole.2.gain <- rownames(sample.information.80x)[sample.information.80x$`2.gain`=="T" & sample.information.80x$`2p.gain`=="F"]
tumors.with.11q.loss <- rownames(sample.information.80x)[sample.information.80x$`11q.deletion`=="T"]

mat <- rbind(mat, matrix("", nrow=10, ncol=ncol(mat),
                         dimnames=list(c("17q gain", "17 gain", "7q gain", "7 gain", "1q gain", "1 gain", "1p loss",
                                       "2 gain", "2p gain", "11q loss"), colnames(mat))))


mat["17q gain", tumors.with.17q] <- "Gain"
mat["17 gain", tumors.with.whole.17.gain] <- "Gain"
mat["7q gain", tumors.with.7q] <- "Gain"
mat["7 gain", tumors.with.whole.7.gain] <- "Gain"
mat["1q gain", tumors.with.1q.gain] <- "Gain"
mat["1 gain", tumors.with.whole.1.gain] <- "Gain"
mat["2p gain", tumors.with.2p] <- "Gain"
mat["2 gain", tumors.with.whole.2.gain] <- "Gain"
mat["17q gain", tumors.with.17q] <- "Gain"
mat["11q loss", tumors.with.11q.loss] <- "Loss"
mat["1p loss", tumors.with.1p.loss] <- "Loss"



mat[is.na(mat)] <- ""
mat[mat %in% c(";wt", ";NA")] <- ""
mat[mat=="whole chromosome loss; stopgain"] <- "stopgain"

## show whole chromosome info only if no partial gain
mat["17 gain", mat["17q gain",]!=""] <- ""
mat["1 gain", mat["1q gain",]!=""] <- ""
mat["2 gain", mat["2p gain",]!=""] <- ""
mat["7 gain", mat["7q gain",]!=""] <- ""