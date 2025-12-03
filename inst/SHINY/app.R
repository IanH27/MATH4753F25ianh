# app.R
# Shiny MLE Explorer: Normal, Exponential, Poisson, Bernoulli, Gamma, Lognormal
# Author: ChatGPT

if (!requireNamespace("shiny", quietly = TRUE)) stop("Please install.packages('shiny')")
library(shiny)

# ---------- Utility: MLE functions ----------
mle_normal <- function(x) {
  mu_hat <- mean(x)
  sigma_hat <- sqrt(mean((x - mu_hat)^2))  # MLE uses 1/n, not 1/(n-1)
  list(par = c(mu = mu_hat, sigma = sigma_hat))
}

mle_exponential <- function(x) {
  # parameterization: rate lambda
  lambda_hat <- 1/mean(x)
  list(par = c(lambda = lambda_hat))
}

mle_poisson <- function(x) {
  lambda_hat <- mean(x)
  list(par = c(lambda = lambda_hat))
}

mle_bernoulli <- function(x) {
  # x in {0,1}
  p_hat <- mean(x)
  # constrain into (0,1) for stability
  p_hat <- min(max(p_hat, 1e-8), 1 - 1e-8)
  list(par = c(p = p_hat))
}

mle_lognormal <- function(x) {
  lx <- log(x)
  mu_hat <- mean(lx)
  sigma_hat <- sqrt(mean((lx - mu_hat)^2))
  list(par = c(meanlog = mu_hat, sdlog = sigma_hat))
}

mle_gamma <- function(x) {
  # Parameterization: shape (k) and rate
  # Use log-params to enforce positivity in optim
  nll <- function(tpars) {
    shape <- exp(tpars[1]); rate <- exp(tpars[2])
    # negative log-likelihood
    -(sum(dgamma(x, shape = shape, rate = rate, log = TRUE)))
  }
  # good starts from method-of-moments
  m <- mean(x); v <- var(x)
  shape0 <- max(m^2 / v, 1e-3)
  rate0  <- max(shape0 / m, 1e-3)
  fit <- optim(par = log(c(shape0, rate0)), fn = nll, method = "BFGS", hessian = TRUE)
  est <- exp(fit$par)
  names(est) <- c("shape", "rate")
  list(par = est, optim = fit)
}

# Distribution metadata for UI and plotting
dist_defs <- list(
  "Normal" = list(
    kind = "continuous",
    params = list(mu = 0, sigma = 1),
    ui = function(inputId) tagList(
      numericInput(paste0(inputId,"_mu"), "True mean (μ)", value = 0),
      sliderInput(paste0(inputId,"_sigma"), "True sd (σ)", min = 0.05, max = 5, value = 1, step = 0.05)
    ),
    rfun = function(n, par) rnorm(n, mean = par["mu"], sd = par["sigma"]),
    d_true = function(x, par) dnorm(x, mean = par["mu"], sd = par["sigma"]),
    d_mle  = function(x, theta) dnorm(x, mean = theta["mu"], sd = theta["sigma"]),
    mle = mle_normal,
    n_params = 2
  ),
  "Exponential" = list(
    kind = "continuous",
    params = list(lambda = 1),
    ui = function(inputId) tagList(
      sliderInput(paste0(inputId,"_lambda"), "True rate (λ)", min = 0.05, max = 5, value = 1, step = 0.05)
    ),
    rfun = function(n, par) rexp(n, rate = par["lambda"]),
    d_true = function(x, par) dexp(x, rate = par["lambda"]),
    d_mle  = function(x, theta) dexp(x, rate = theta["lambda"]),
    mle = mle_exponential,
    n_params = 1
  ),
  "Poisson" = list(
    kind = "discrete",
    params = list(lambda = 3),
    ui = function(inputId) tagList(
      sliderInput(paste0(inputId,"_lambda"), "True mean (λ)", min = 0.2, max = 15, value = 3, step = 0.2)
    ),
    rfun = function(n, par) rpois(n, lambda = par["lambda"]),
    pmf_true = function(k, par) dpois(k, lambda = par["lambda"]),
    pmf_mle  = function(k, theta) dpois(k, lambda = theta["lambda"]),
    mle = mle_poisson,
    n_params = 1
  ),
  "Bernoulli" = list(
    kind = "discrete",
    params = list(p = 0.5),
    ui = function(inputId) tagList(
      sliderInput(paste0(inputId,"_p"), "True success prob (p)", min = 0.01, max = 0.99, value = 0.5, step = 0.01)
    ),
    rfun = function(n, par) rbinom(n, size = 1, prob = par["p"]),
    pmf_true = function(k, par) dbinom(k, size = 1, prob = par["p"]),
    pmf_mle  = function(k, theta) dbinom(k, size = 1, prob = theta["p"]),
    mle = mle_bernoulli,
    n_params = 1
  ),
  "Gamma" = list(
    kind = "continuous",
    params = list(shape = 2, rate = 1),
    ui = function(inputId) tagList(
      sliderInput(paste0(inputId,"_shape"), "True shape", min = 0.3, max = 8, value = 2, step = 0.1),
      sliderInput(paste0(inputId,"_rate"),  "True rate",  min = 0.1, max = 5, value = 1, step = 0.05)
    ),
    rfun = function(n, par) rgamma(n, shape = par["shape"], rate = par["rate"]),
    d_true = function(x, par) dgamma(x, shape = par["shape"], rate = par["rate"]),
    d_mle  = function(x, theta) dgamma(x, shape = theta["shape"], rate = theta["rate"]),
    mle = mle_gamma,
    n_params = 2
  ),
  "Lognormal" = list(
    kind = "continuous",
    params = list(meanlog = 0, sdlog = 0.6),
    ui = function(inputId) tagList(
      numericInput(paste0(inputId,"_meanlog"), "True meanlog", value = 0),
      sliderInput(paste0(inputId,"_sdlog"), "True sdlog", min = 0.05, max = 2, value = 0.6, step = 0.05)
    ),
    rfun = function(n, par) rlnorm(n, meanlog = par["meanlog"], sdlog = par["sdlog"]),
    d_true = function(x, par) dlnorm(x, meanlog = par["meanlog"], sdlog = par["sdlog"]),
    d_mle  = function(x, theta) dlnorm(x, meanlog = theta["meanlog"], sdlog = theta["sdlog"]),
    mle = mle_lognormal,
    n_params = 2
  )
)

# ---------- UI ----------
ui <- fluidPage(
  titlePanel("Maximum Likelihood Estimation Demo (Univariate Distributions)"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dist", "Distribution",
                  choices = names(dist_defs), selected = "Normal"),
      uiOutput("param_ui"),
      hr(),
      sliderInput("n", "Sample size (n)", min = 5, max = 2000, value = 200, step = 5),
      numericInput("seed", "Random seed", value = 123, min = 1),
      actionButton("resample", "Generate / Refit", class = "btn-primary"),
      hr(),
      helpText("Tip: Adjust parameters and sample size, then click Generate / Refit to see MLE in action.")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Data & Fits",
                 plotOutput("fit_plot", height = 420),
                 tableOutput("mle_table")
        ),
        tabPanel("Log-Likelihood",
                 conditionalPanel(
                   "['Exponential','Poisson','Bernoulli'].indexOf(input.dist) >= 0",
                   plotOutput("ll_profile", height = 420),
                   helpText("Log-likelihood profile over the single parameter.")
                 ),
                 conditionalPanel(
                   "['Normal','Gamma','Lognormal'].indexOf(input.dist) >= 0",
                   plotOutput("ll_contour", height = 480),
                   helpText("Log-likelihood contour over two parameters (centered around the MLE).")
                 )
        ),
        tabPanel("About",
                 tags$h4("What you're seeing"),
                 tags$p("This app simulates data from a chosen distribution with user-selected 'true' parameters, ",
                        "then computes the maximum likelihood estimates (MLEs). For continuous distributions, we overlay ",
                        "the true density and the MLE-fitted density on a histogram. For discrete distributions, we show a bar plot ",
                        "of empirical probabilities, with true and MLE probability mass overlays."),
                 tags$ul(
                   tags$li("Closed-form MLEs: Normal (μ, σ), Exponential (λ), Poisson (λ), Bernoulli (p), Lognormal (meanlog, sdlog)."),
                   tags$li("Numerical MLE (via optim): Gamma (shape, rate).")
                 ),
                 tags$p("The log-likelihood tab illustrates how the MLE maximizes the log-likelihood: ",
                        "a 1D profile for 1-parameter models and a 2D contour for 2-parameter models.")
        )
      )
    )
  )
)

# ---------- SERVER ----------
server <- function(input, output, session) {
  
  # Dynamic parameter UI
  output$param_ui <- renderUI({
    dist_defs[[input$dist]]$ui("par")
  })
  
  # Gather "true" parameter vector from UI inputs
  true_params <- reactive({
    d <- input$dist
    def <- dist_defs[[d]]
    nm <- names(def$params)
    vals <- setNames(numeric(length(nm)), nm)
    for (k in nm) {
      vals[k] <- input[[paste0("par_", k)]]
    }
    vals
  })
  
  # Data generation, re-run when clicking button or changing seed/dist/params
  data_r <- eventReactive(input$resample, {
    set.seed(input$seed)
    d <- input$dist
    par <- true_params()
    n <- input$n
    x <- dist_defs[[d]]$rfun(n, par)
    list(x = x, par = par, dist = d)
  }, ignoreInit = TRUE)
  
  # Compute MLE
  mle_r <- reactive({
    req(data_r())
    x <- data_r()$x
    d <- data_r()$dist
    fit <- dist_defs[[d]]$mle(x)
    fit$par
  })
  
  # --- Plot: Data with true vs MLE overlays ---
  output$fit_plot <- renderPlot({
    req(data_r())
    x <- data_r()$x
    d <- data_r()$dist
    par_true <- data_r()$par
    par_hat  <- mle_r()
    
    if (dist_defs[[d]]$kind == "continuous") {
      # histogram with density
      hist(x, breaks = "FD", freq = FALSE, main = paste0(d, ": Data & Fits"),
           xlab = "x", border = "gray50", col = "white")
      rng <- range(x)
      xs <- seq(rng[1] - 0.1*diff(rng), rng[2] + 0.1*diff(rng), length.out = 400)
      
      # true density
      lines(xs, dist_defs[[d]]$d_true(xs, par_true), lwd = 2, lty = 1)
      # MLE density
      lines(xs, dist_defs[[d]]$d_mle(xs, par_hat), lwd = 2, lty = 2)
      
      legend("topright", bty = "n", lwd = 2, lty = c(1,2),
             legend = c("True density", "MLE density"))
    } else {
      # discrete: bar plot of empirical probabilities
      tab <- table(x)
      ks <- as.integer(names(tab))
      probs_emp <- as.numeric(tab) / length(x)
      plot(ks, probs_emp, type = "h", lwd = 6, lend = 2,
           xlab = "k", ylab = "Probability", main = paste0(d, ": Empirical vs True/MLE PMF"),
           ylim = c(0, max(probs_emp)*1.2))
      
      # True PMF
      k_grid <- seq(min(ks), max(ks), by = 1)
      points(k_grid, dist_defs[[d]]$pmf_true(k_grid, par_true), type = "b", pch = 19)
      
      # MLE PMF
      points(k_grid, dist_defs[[d]]$pmf_mle(k_grid, par_hat), type = "b", pch = 1, lty = 2)
      
      legend("topright", bty = "n",
             legend = c("Empirical", "True PMF", "MLE PMF"),
             pch = c(NA, 19, 1), lty = c(1,1,2), lwd = c(6,1,1),
             pt.cex = c(NA, 1, 1))
    }
  })
  
  # --- MLE summary table ---
  output$mle_table <- renderTable({
    req(data_r())
    d <- data_r()$dist
    par_true <- data_r()$par
    par_hat <- mle_r()
    
    # closed-form notes
    notes <- switch(d,
                    "Normal"     = "Closed-form: μ̂=mean(x), σ̂=sqrt(mean((x-μ̂)^2))",
                    "Exponential"= "Closed-form: λ̂=1/mean(x)",
                    "Poisson"    = "Closed-form: λ̂=mean(x)",
                    "Bernoulli"  = "Closed-form: p̂=mean(x)",
                    "Lognormal"  = "Closed-form: meanloĝ=mean(log x), sdloĝ=sd(log x) with 1/n",
                    "Gamma"      = "Numerical MLE via optim (shape, rate)"
    )
    
    data.frame(
      Parameter = union(names(par_true), names(par_hat)),
      True = signif(par_true[match(union(names(par_true), names(par_hat)), names(par_true))], 6),
      MLE  = signif(par_hat[match(union(names(par_true), names(par_hat)), names(par_hat))], 6),
      Note = c(notes, rep("", length(union(names(par_true), names(par_hat))) - 1)),
      check.names = FALSE
    )
  })
  
  # ---------- Log-likelihood (1-parameter) ----------
  output$ll_profile <- renderPlot({
    req(data_r())
    x <- data_r()$x
    d <- data_r()$dist
    par_hat <- mle_r()
    if (!(d %in% c("Exponential","Poisson","Bernoulli"))) return(invisible())
    
    if (d == "Exponential") {
      # lambda in (0, ∞)
      lam_hat <- par_hat["lambda"]
      grid <- seq(max(1e-4, lam_hat/5), lam_hat*5, length.out = 300)
      ll <- sapply(grid, function(l) sum(dexp(x, rate = l, log = TRUE)))
      plot(grid, ll, type = "l", lwd = 2, xlab = "λ", ylab = "log-likelihood",
           main = "Exponential: log-likelihood vs λ")
      abline(v = lam_hat, lty = 2)
    } else if (d == "Poisson") {
      lam_hat <- par_hat["lambda"]
      grid <- seq(max(1e-4, lam_hat/5), lam_hat*5, length.out = 300)
      ll <- sapply(grid, function(l) sum(dpois(x, lambda = l, log = TRUE)))
      plot(grid, ll, type = "l", lwd = 2, xlab = "λ", ylab = "log-likelihood",
           main = "Poisson: log-likelihood vs λ")
      abline(v = lam_hat, lty = 2)
    } else if (d == "Bernoulli") {
      p_hat <- par_hat["p"]
      grid <- seq(1e-4, 1-1e-4, length.out = 300)
      ll <- sapply(grid, function(p) sum(dbinom(x, size = 1, prob = p, log = TRUE)))
      plot(grid, ll, type = "l", lwd = 2, xlab = "p", ylab = "log-likelihood",
           main = "Bernoulli: log-likelihood vs p")
      abline(v = p_hat, lty = 2)
    }
  })
  
  # ---------- Log-likelihood contour (2-parameter) ----------
  output$ll_contour <- renderPlot({
    req(data_r())
    x <- data_r()$x
    d <- data_r()$dist
    par_hat <- mle_r()
    if (!(d %in% c("Normal","Gamma","Lognormal"))) return(invisible())
    
    if (d == "Normal") {
      mu_hat <- par_hat["mu"]; sigma_hat <- par_hat["sigma"]
      mu_grid <- seq(mu_hat - 2*sd(x), mu_hat + 2*sd(x), length.out = 60)
      sigma_grid <- seq(max(0.05, sigma_hat/3), sigma_hat*3, length.out = 60)
      z <- outer(mu_grid, sigma_grid, Vectorize(function(m, s) sum(dnorm(x, m, s, log = TRUE))))
      contour(mu_grid, sigma_grid, z,
              xlab = "μ", ylab = "σ", main = "Normal: log-likelihood contour")
      points(mu_hat, sigma_hat, pch = 19)
      legend("topright", bty = "n", pch = 19, legend = "MLE")
    }
    
    if (d == "Gamma") {
      shape_hat <- par_hat["shape"]; rate_hat <- par_hat["rate"]
      shp_grid <- seq(max(0.2, shape_hat/3), shape_hat*3, length.out = 60)
      rate_grid <- seq(max(0.05, rate_hat/3), rate_hat*3, length.out = 60)
      z <- outer(shp_grid, rate_grid, Vectorize(function(a, r) sum(dgamma(x, shape = a, rate = r, log = TRUE))))
      contour(shp_grid, rate_grid, z,
              xlab = "shape", ylab = "rate", main = "Gamma: log-likelihood contour")
      points(shape_hat, rate_hat, pch = 19)
      legend("topright", bty = "n", pch = 19, legend = "MLE")
    }
    
    if (d == "Lognormal") {
      ml_hat <- par_hat["meanlog"]; sl_hat <- par_hat["sdlog"]
      ml_grid <- seq(ml_hat - 2*sl_hat, ml_hat + 2*sl_hat, length.out = 60)
      sl_grid <- seq(max(0.05, sl_hat/3), sl_hat*3, length.out = 60)
      z <- outer(ml_grid, sl_grid, Vectorize(function(m, s) sum(dlnorm(x, meanlog = m, sdlog = s, log = TRUE))))
      contour(ml_grid, sl_grid, z,
              xlab = "meanlog", ylab = "sdlog", main = "Lognormal: log-likelihood contour")
      points(ml_hat, sl_hat, pch = 19)
      legend("topright", bty = "n", pch = 19, legend = "MLE")
    }
  })
}

shinyApp(ui, server)