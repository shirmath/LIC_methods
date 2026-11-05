# H3: both methods are within 20% of sigma on a smooth signal at 48 frames; diff2's bias on 6 cycles is pinned

    Code
      cat("diff2 bias (estimate / sigma) on a 6-cycle sinusoid at noise ratio 0.5:",
        round(mean(est_fast) / sigma, 2))
    Output
      diff2 bias (estimate / sigma) on a 6-cycle sinusoid at noise ratio 0.5: 1.11

# R5: supports along the path are nearly nested [snapshot of the non-nested count]

    Code
      cat("non-nested support steps over a 30-point path (seed 101, n_D = 36):",
        non_nested)
    Output
      non-nested support steps over a 30-point path (seed 101, n_D = 36): 2

