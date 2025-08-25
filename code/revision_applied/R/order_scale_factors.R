

order_and_scale <- function(F_hat){
    # sort factors
    i <- apply(abs(F_hat), 2, which.max)
    o <- order(i)
    F_hat <- F_hat[, o]
    # flip so that largest element is positive
    i <- i[o]
    fct_sign <- sapply(1:length(i), function(x){sign(F_hat[i[x], x])})
    F_hat <- t(t(F_hat)*fct_sign)
    return(list("F_hat" = F_hat, order = o))
}



