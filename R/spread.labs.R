# Copied from the orphaned package TeachingDemos version 2.12.1 on Feb 16, 2024.
# ------------------------------------------------------------------------------
#
# --Title--
#
# Spread out close points for labeling in plots
#
# --Description--
#
# This function takes as set of coordinates and spreads out the close
# values so that they can be used in labeling plots without overlapping.
#
# --Usage--
#
# spread.labs(x, mindiff, maxiter = 1000, stepsize = 1/10, min = -Inf, max = Inf)
#
# --Arguments--
#
#   x         The coordinate values (x or y, not both) to spread out.
#   mindiff   The minimum distance between return values
#   maxiter   The maximum number of iterations
#   stepsize  How far to move values in each iteration
#   min       Minimum bound for returned values
#   max       Maximum bound for returned values
#
# --Details--
#
#   Sometimes the desired locations for labels in plots results in the
#   labels overlapping.  This function takes the coordinate values (x or
#-  y, not both) and finds those points that are less than mindiff
#   (usually a function of strheight or strwidth ) apart and
#   increases the space between them (by stepsize * mindiff ).
#   This may or may not be enough and moving some points
#   away from their nearest neighbor may move them too close to another
#   neighbor, so the process is iterated until either maxiter steps
#   have been tried, or all the values are at least mindiff apart.
#
#   The min and max arguments prevent the values from going
#   outside that range (they should be specified such that the original
#   values are all inside the range).
#
#   The values do not need to be presorted.
#
# --Return Value--
#
#   A vector of coordinates (order corresponding to the original x )
#   that can be used as a replacement for x in placing labels.
#
# --Author--
#
# Greg Snow  email 538280@gmail.com
#
# --See Also--
#
# The spread.labels function in the plotrix package.
#
# --Examples--
#
# # overlapping labels
# plot(as.integer(state.region), state.x77[,1], ylab='Population',
#     xlab='Region',xlim=c(1,4.75), xaxt='n')
# axis(1, at=1:4, lab=levels(state.region) )
#
# text( as.integer(state.region)+.5, state.x77[,1], state.abb )
# segments( as.integer(state.region)+0.025, state.x77[,1],
#         as.integer(state.region)+.375, state.x77[,1] )
#
# # now lets redo the plot without overlap
#
# tmp.y <- state.x77[,1]
# for(i in levels(state.region) ) {
#   tmp <- state.region == i
#   tmp.y[ tmp ] <- spread.labs( tmp.y[ tmp ], 1.2*strheight('A'),
#     maxiter=1000, min=0 )
# }
#
# plot(as.integer(state.region), state.x77[,1], ylab='Population',
#     xlab='Region', xlim=c(1,4.75), xaxt='n')
# axis(1, at=1:4, lab=levels(state.region) )
#
# text( as.integer(state.region)+0.5, tmp.y, state.abb )
# segments( as.integer(state.region)+0.025, state.x77[,1],
#         as.integer(state.region)+0.375, tmp.y )
# }

spread.labs <- function(x, mindiff, maxiter=1000, stepsize=1/10,
                          min=-Inf, max=Inf) {
    unsort <- order(order(x))
    x <- sort(x)
    df <- x[-1] - x[ -length(x) ]

    stp <- mindiff * stepsize

    i <- 1
    while( any( df < mindiff ) ) {
        tmp <- c( df < mindiff, FALSE )
        if( tmp[1] && (x[1] - stp) < min ) {  # don't move bottom set
            tmp2 <- as.logical( cumprod(tmp) )
            tmp <- tmp & !tmp2
        }
        x[ tmp ] <- x[ tmp ] - stp
        tmp <- c( FALSE, df < mindiff )
        if( tmp[length(tmp)] && (x[length(x)] + stp) > max ) { # don't move top
            tmp2 <- rev( as.logical( cumprod( rev(tmp) ) ) )
            tmp <- tmp & !tmp2
        }
        x[ tmp ] <- x[ tmp] + stp

        df <- x[-1] - x[-length(x)]
        i <- i + 1
        if( i > maxiter ) {
            warning("Maximum iterations reached")
            break
        }
    }
    x[unsort]
}
