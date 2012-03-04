\name{cqnnews}
\title{CQN News}
\encoding{UTF-8}

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
