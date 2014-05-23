getPValue <- function( ## calculate p-value from an element of type 'history'
                      history,
                      nObs){
  y <- history[nrow(history), ]
  psi <- y[["psi"]]
  phi <- y[["phi"]]
  se <- y[["psi.sd"]]/sqrt(nObs)
  1-pnorm(abs(psi-phi), sd=se)
}
