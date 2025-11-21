spaceout_points <- function(
        x,
        limit = NULL,
        margin = 0.05,
        spacing = NULL,
        iterations = 500,
        force = 1,
        force_decay = 0.0,
        attraction = 0.05,
        attr_decay = 0.0,
        ret_history = F) {
    ixs <- order(x)
    x <- x[ixs]
    n <- length(x)
    r <- range(x)
    s <- diff(r)
    if (is.null(limit)) {
        start <- r[1] - margin * s
        end <- r[2] + margin * s
    } else {
        start <- limit[1]
        end <- limit[2]
        stopifnot(start < r[1] & r[2] < end)
    }
    if (is.null(spacing)) {
        spacing <- (end - start) / (n - 1)
    }
    eps <- s * 1e-10
    
    if (ret_history) {
        x_history <- data.frame(x = x, transformed = x, iter = '0')
    }
    
    curr_x <- x
    curr_force <- force
    curr_attr <- attraction
    for (i in 1:iterations) {
        x_ext <- c(start, curr_x, end)
        d <- diff(x_ext)
        force_val <- spacing * curr_force * exp(-d / spacing)
        sum_force <- force_val[1:n] - force_val[2:(n + 1)]
        x_next1 <- curr_x + sum_force
        
        d_init <- x_next1 - x
        attr_val <- pmin(curr_attr * d_init, spacing) |> pmax(-spacing)
        x_next2 <- x_next1 - attr_val

        midpoints <- 0.5 * (curr_x[1:(n-1)] + curr_x[2:n])
        midpoints_l <- c(start, midpoints + eps)
        midpoints_r <- c(midpoints - eps, end)
        x_next <- x_next2 |>
            pmin(midpoints_r - eps) |>
            pmax(midpoints_l + eps)

        if (ret_history) {
            x_history <- rbind(x_history,
                data.frame(x, transformed = x_next, iter = as.character(i)))
        }

        d <- x_next - curr_x
        curr_x <- x_next
        if (max(abs(d)) < eps) {
            break
        }
        curr_force <- curr_force * (1 - force_decay)
        curr_attr <- curr_attr * (1 - attr_decay)
    }
    curr_x[1] <- pmin(curr_x[1], x[1])
    curr_x[n] <- pmax(curr_x[n], x[n])
    if (ret_history) {
        rbind(x_history, data.frame(x = x, transformed = curr_x, iter = 'final')) |>
            dplyr::mutate(iter = factor(iter, levels = unique(iter)))
    } else {
        curr_x[order(ixs)]
    }
}
