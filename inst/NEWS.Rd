\name{cqnnews}
\title{CQN News}
\encoding{UTF-8}

\section{Version 1.3}{
  \itemize{
    \item Bugfix to the vignette; the two-panel (color) plot on page 6
    used CQN corrected data as blue points in both panels.  Now the left
    plot shows standard RPKM in blue.  Thanks to Maria Keays
    \email{mkeays@ebi.ac.uk}. 
    \item Small fix to the vignette in the edgeR example, caused by
    changes to edgeR.
    \item Updated the citations in the vignette and the CITATION file.
    \item cqn.fixedlength has been removed, using cqn(lengthMethod =
    "fixed") instead.
    \item A call slot has been added to cqn objects.
    \item cqn() now accepts count matrices with 1 column or vectors (although it makes
    little sense to use the function on such data).
    \item Added Questions and Answers to the vignette and moved vignette
    to vignettes dir.
  }
}


\section{Version 1.1}{
  \itemize{
    \item Same fixes as in 1.0.1
    \item Resaved the data files, so they take up less space.
    \item Added edgeR as a Suggests: since the vignette uses it.
  }
}

\section{Version 1.0.1}{
  \itemize{
    \item The function alpha has been moved from ggplot2 to the new
      package scales.  Vignette and Suggests: fields have been changed
      accordingly.
    \item sizeFactors = NULL should now work for cqn and cq.fixedlength.  
      Thanks to Maria Chikina for reporting this
  }
}

\section{Version 0.99}{
  \itemize{
    \item Version bump as well as added to Bioconductor devel branch.
    \item Wrote vignette.
    \item Added additional manpages as well as examples.
    \item Added data from Montgomery et al.
    \item Added NAMESPACE file.
    \item Added manpage for cqn.
    \item Renamed cqn, cqn2.
    \item Added NEWS.Rd file.
    \item Added vignette skeleton.
    \item Added CITATION file.
  }
}
