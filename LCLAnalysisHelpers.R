# LCLAnalysisHelpers.R
# Copyright (C) 2020 Lestat Ali
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

library(metaseqR)

apply.fisher.method <- function (pvalues1, pvalues2, zero.sub = 1e-05) {
  pvals <- matrix(c(pvalues1, pvalues2), ncol = 2)
  pvals[pvals == 0] <- zero.sub
  fisher.sums <- data.frame(do.call(rbind, apply(pvals, 1, fisher.sum, zero.sub = zero.sub, na.rm = F)))
  rownames(fisher.sums) <- rownames(pvals)
  fisher.sums$p.value <- pchisq(fisher.sums$S, df = 2 * fisher.sums$num.p, lower.tail = F)
  fisher.sums$p.adj <- p.adjust(fisher.sums$p.value, "BH")
  return(fisher.sums)
}
