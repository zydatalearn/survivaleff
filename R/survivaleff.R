####变量基本特征####
#默认保留两位小数
#' According to the provided single variable, automatically identify whether it is a categorical variable
#' or a continuous variable. If it is a continuous variable,
#' provide the basic characteristics of this variable and conduct a normality test.
#' If it is a categorical variable, provide the frequency and percentage.
#'
#'
#' @title summary_variable
#' @param data data object.
#' @param variable single variable.
#' @author camellia
#' @examples
#' set.seed(123)
#' dt <- data.frame(
#' age = sample(20:80,100,replace = TRUE),
#' group = sample(c("A", "B"),100,replace = TRUE))
#' summary_variable(data = dt, col = "age", digits = 2,p_digits=3)
#' summary_variable(data = dt, col = "group", digits = 2,p_digits=3)
#' @return A part of input data.
#' @export
summary_variable <- function(data, cols, digits = 2, p_digits = 3) {
  options(scipen = 999)

  format_p <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 10^(-p_digits)) {
      return(paste0("p < ", formatC(10^(-p_digits), format = "f", digits = p_digits)))
    } else {
      return(paste0("p = ", formatC(p, format = "f", digits = p_digits)))
    }
  }

  all_names <- names(data)
  cols_resolved <- sapply(cols, function(col) {
    if (is.numeric(col)) {
      if (col > 0 && col <= ncol(data)) return(all_names[col]) else return(NA)
    } else if (is.character(col)) {
      if (col %in% all_names) return(col) else return(NA)
    } else {
      return(NA)
    }
  }, USE.NAMES = FALSE)

  if (any(is.na(cols_resolved))) {
    warning("以下变量在数据中未找到并已跳过: ",
            paste(cols[is.na(cols_resolved)], collapse = ", "))
  }

  valid_cols <- na.omit(cols_resolved)

  results_num <- data.frame(
    变量名 = character(),
    样本数 = numeric(),
    `均值±标准差` = character(),
    `中位数(最小, 最大)` = character(),
    `中位数(P25, P75)` = character(),
    Shapiro_W = character(),
    Shapiro_p = character(),
    stringsAsFactors = FALSE
  )

  results_cat <- data.frame(
    变量名 = character(),
    亚组名称 = character(),
    `频数(百分比)` = character(),
    stringsAsFactors = FALSE
  )

  for (colname in valid_cols) {
    x <- data[[colname]]
    na_count <- sum(is.na(x))
    x_no_na <- na.omit(x)
    n <- length(x_no_na)

    if (is.numeric(x)) {
      quartiles <- quantile(x_no_na, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
      mean_val <- mean(x_no_na)
      sd_val <- sd(x_no_na)
      median_val <- median(x_no_na)

      mean_sd <- sprintf(paste0("%.", digits, "f ± %.", digits, "f"), mean_val, sd_val)
      median_range <- sprintf(paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"),
                              median_val, quartiles[1], quartiles[5])
      median_iqr <- sprintf(paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"),
                            median_val, quartiles[2], quartiles[4])

      cat("\n--- 数值变量: ", colname, " ---\n")
      cat("样本量（非缺失）:", n, "\n")
      cat("缺失值数量:", na_count, "\n")
      cat(sprintf("最小值: %.2f  第一四分位数: %.2f  中位数: %.2f  第三四分位数: %.2f  最大值: %.2f\n",
                  quartiles[1], quartiles[2], quartiles[3], quartiles[4], quartiles[5]))
      cat(sprintf("均值: %.2f  标准差: %.2f\n", mean_val, sd_val))

      if (n >= 3 && n <= 5000) {
        shapiro_res <- shapiro.test(x_no_na)
        W <- round(shapiro_res$statistic, digits)
        p_val <- shapiro_res$p.value
        cat("Shapiro-Wilk 检验结果:\n")
        cat("W =", W, " | ", format_p(p_val), "\n")
        print(shapiro_res)
      } else {
        W <- NA
        p_val <- NA
        cat("样本量不足，跳过正态性检验\n")
      }

      results_num <- rbind(
        results_num,
        data.frame(
          变量名 = colname,
          样本数 = n,
          `均值±标准差` = mean_sd,
          `中位数(最小, 最大)` = median_range,
          `中位数(P25, P75)` = median_iqr,
          Shapiro_W = ifelse(is.na(W), "NA", sprintf(paste0("%.", digits, "f"), W)),
          Shapiro_p = format_p(p_val),
          stringsAsFactors = FALSE
        )
      )
    } else {
      cat("\n--- 分类变量: ", colname, " ---\n")
      cat("缺失值数量:", na_count, "\n")
      tab <- table(x_no_na, useNA = "no")
      tab_prop <- prop.table(tab)

      for (lvl in names(tab)) {
        cnt <- tab[[lvl]]
        pct <- sprintf(paste0("%.", digits, "f%%"), tab_prop[[lvl]] * 100)
        combined <- paste0(cnt, " (", pct, ")")
        cat(sprintf("%s: %s\n", lvl, combined))

        results_cat <- rbind(
          results_cat,
          data.frame(
            变量名 = colname,
            亚组名称 = lvl,
            `频数(百分比)` = combined,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }

  if (nrow(results_num) > 0 && nrow(results_cat) > 0) {
    return(list(数值变量 = results_num, 分类变量 = results_cat))
  } else if (nrow(results_num) > 0) {
    return(results_num)
  } else if (nrow(results_cat) > 0) {
    return(results_cat)
  } else {
    warning("未找到符合条件的变量")
    return(NULL)
  }
}
####生存曲线基本特征-ALL####
#' When drawing the overall survival curve, the median time and 95% CI confidence
#' interval can be quickly viewed and the survival rates
#' at the specified times can be output.
#'
#'
#' @title survival_overview
#' @param data data object.
#' @param time_var time variable.
#' @param status_var Numeric type, only 0/1,0 is censored and 1 is the occurrence of an event.
#' @param time_points Specific 3 time points.
#' @param digits Retain the number of decimal places, default is 2
#' @author camellia
#' @examples
#' library(survival)
#' set.seed(123)
#' n <- 100
#' dt <- data.frame(
#' pfs1 = rexp(n, rate = 0.1),
#' pfs1_status = sample(c(0, 1), n, replace = TRUE))
#' survival_overview(data = dt, time_var = "pfs1", status_var = "pfs1_status", time_points = c(6, 12, 24))
#' @return A part of input data.
#' @export

survival_overview <- function(data, time_var, status_var,time_points=c(6,12),digits=2) {
  time_var <- data[[time_var]]
  status_var <- data[[status_var]]
  fit <- survfit(Surv(time_var,status_var)~1, data = data)
  valid_idx <- !is.na(time_var) & !is.na(status_var)
  time_var <- time_var[valid_idx]
  total_count <- sum(valid_idx)
  event_count <- sum(status_var == 1,na.rm = T)
  event_percentage <- sprintf(paste0('%.', digits, 'f'),event_count / total_count * 100)

  event_info <- paste0("Overall: ", event_count, "/", total_count, " (", event_percentage, "%)")

  mpfs <- paste0(sprintf(paste0('%.', digits, 'f'), quantile(fit)$quantile[[2]]),
                 " month", "(", sprintf(paste0('%.', digits, 'f'), quantile(fit)$lower[[2]]),
                 "-", sprintf(paste0('%.', digits, 'f'), quantile(fit)$upper[[2]]), ")")
  surv_values <- lapply(time_points, function(t) {
    surv <- sprintf(paste0('%.', digits, 'f'), summary(fit, times = c(t))$surv[1] * 100)
    lower <- sprintf(paste0('%.', digits, 'f'), summary(fit, times = c(t))$lower[1]* 100)
    upper <- sprintf(paste0('%.', digits, 'f'), summary(fit, times = c(t))$upper[1]* 100)
    surv_text <- paste0(t,"months:",surv,"(95%CI:",lower,"-",upper,")")
    return(surv_text)
  })
  names(surv_values) <- paste0(time_points, " months")
  result <- list(
    KM_imformation = fit,
    Event_rate = event_info,
    Median_time = mpfs,
    surv_rate = surv_values
  )
  return(result)
  cat("\nKM information:\n")
  print(fit)
  cat("\nEvent/Total:\n")
  print(event_info)
  cat("\nMedian(95%CI):\n")
  print(mpfs)
  cat("\nsurvival rate/month:\n")
  print(surv_values)
}

####亚组生存曲线基本信息:HR,HR_pvalue,Log-rank,Median####
#' When conducting subgroup analysis, the model characteristics can be quickly obtained,
#' including the median survival data of each group,
#' subgroup results, HR value, and Log-rank test results.
#'
#'
#' @title survival_analysis_by_group
#' @param data data object.
#' @param time_var time variable.
#' @param status_var Numeric type, only 0/1,0 is censored and 1 is the occurrence of an event.
#' @param group_var subgroup variable..
#' @author camellia
#' @examples
#' library(survival)
#' set.seed(123)
#' dt <- data.frame(
#' pfs1 = c(5, 10, 15, 20, 30, 35, 40, 50, 60, 70),
#' pfs1_status = c(1, 1, 0, 1, 0, 1, 1, 0, 1, 1),
#' 肝转移 = factor(c("Yes", "No", "Yes", "No", "Yes", "No", "Yes", "No", "Yes", "No")))
#' result <- survival_analysis_by_group(dt, "pfs1", "pfs1_status", "肝转移")
#' result
#' @return A part of input data.
#' @export
survival_analysis_by_group <- function(data, time_var, status_var, group_var, digits = 2, p_digits = 3) {
  data <- data.frame(data)
  data[[group_var]] <- factor(data[[group_var]])

  cox_model <- coxph(Surv(data[[time_var]], data[[status_var]]) ~ data[[group_var]], data = data)
  cox_summary <- summary(cox_model)

  HR <- sprintf(paste0('%.', digits, 'f'), cox_summary$coefficients[, 2])
  Pvalue <- ifelse(cox_summary$coefficients[, 5] < 0.001, "<0.001",
                   sprintf(paste0('%.', p_digits, 'f'), cox_summary$coefficients[, 5]))

  CI5 <- sprintf(paste0('%.', digits, 'f'), cox_summary$conf.int[, 3])
  CI95 <- sprintf(paste0('%.', digits, 'f'), cox_summary$conf.int[, 4])
  CI <- paste0(HR, "[95%CI,", CI5, "-", CI95, "]")

  variable_diff <- survdiff(Surv(data[[time_var]], data[[status_var]]) ~ data[[group_var]], data = data)
  p_val <- sprintf(paste0('%.', p_digits, 'f'),1 - pchisq(variable_diff$chisq, length(variable_diff$n) - 1))
  p_val <- ifelse(p_val < 0.001, "<0.001", p_val)

  fit <- survfit(Surv(data[[time_var]], data[[status_var]]) ~ data[[group_var]], data = data)
  PH_test = cox.zph(cox_model)
  fit_summary <- summary(fit)

  group_levels <- levels(data[[group_var]])
  median_survival_by_group <- sapply(group_levels, function(group_name) {
    group_index <- which(group_levels == group_name)
    median_time <- fit_summary$table[group_index, "median"]
    ci_lower <- fit_summary$table[group_index, "0.95LCL"]
    ci_upper <- fit_summary$table[group_index, "0.95UCL"]

    group_data <- data[data[[group_var]] == group_name, ]
    event_count <- sum(group_data[[status_var]] == 1, na.rm = TRUE)
    total_count <- sum(!is.na(group_data[[time_var]]))
    event_percentage <- sprintf(paste0('%.', digits, 'f'),event_count/total_count*100)

    median_time_formatted <- sprintf(paste0('%.', digits, 'f'), median_time)
    ci_lower_formatted <- sprintf(paste0('%.', digits, 'f'), ci_lower)
    ci_upper_formatted <- sprintf(paste0('%.', digits, 'f'), ci_upper)
    event_count_formatted <- event_count
    total_count_formatted <- total_count

    paste0(group_name,"Events:", event_count_formatted, "/",
           total_count_formatted, "(", event_percentage, "%) ",
           ":Median=", median_time_formatted, "(", ci_lower_formatted,
           "-", ci_upper_formatted,
           ")")
  })

  result <- list(
    KM_imformation = fit,
    Cox_Model_Summary = cox_summary,
    PH_test = PH_test,
    HR_CI = CI,
    P_value = Pvalue,
    P_value_Diff = p_val,
    Median_Survival_by_Group = median_survival_by_group
  )
  return(result)
  cat("\nCox model：\n")
  print(cox_summary)
  cat("\nCox model PH test：\n")
  print(PH_test)
  cat("\nKM information：\n")
  print(fit)
  cat("\nthe median of subgroup and 95% CI：\n")
  print(median_survival_by_group)
  cat("\nHR 和 95%CI：\n")
  print(CI)
  cat("\ncox model p-value：\n")
  print(Pvalue)
  cat("\nLog-rank p-value：\n")
  print(p_val)
}


####ORR_DCR计算####
#' Used to quickly calculate ORR, DCR values
#'
#' @title ORR_DCR
#' @param x1 the number of ORR events.
#' @param x2 the number of DCR events.
#' @param n the Total.
#' @author camellia
#' @examples
#' library(DescTools)
#' ORR_DCR(34,40,75)
#' @return A part of input data.
#' @export

ORR_DCR <- function(x1,x2,n,digits=2){
  orr_event_rate <- sprintf(paste0('%.', digits, 'f'),x1/n*100)
  dcr_event_rate <- sprintf(paste0('%.', digits, 'f'),x2/n*100)
  ORR_event_rate = paste0(x1,"/",n,"(",orr_event_rate,"%) ")
  DCR_event_rate = paste0(x2,"/",n,"(",dcr_event_rate,"%) ")
  orr <- BinomCI(x=x1, n=n, method="clopper-pearson")

  orr_rate <- sprintf(paste0('%.', digits, 'f'),orr[[1]]*100)
  orrci <- paste0(sprintf(paste0('%.', digits, 'f'),orr[[2]]*100),"%",
                  "-",
                  sprintf(paste0('%.', digits, 'f'),orr[[3]]*100))

  dcr <- BinomCI(x=x2, n=n, method="clopper-pearson")

  dcr_rate <- sprintf(paste0('%.', digits, 'f'),dcr[[1]]*100)
  dcrci <- paste0(sprintf(paste0('%.', digits, 'f'),dcr[[2]]*100),"%",
                  "-",
                  sprintf(paste0('%.', digits, 'f'),dcr[[3]]*100))

  ORR=paste0(orr_rate,"%(95%CI,",orrci,"%)")
  DCR=paste0(dcr_rate,"%(95%CI,",dcrci,"%)")

  result <- list(
    ORR_event_rate = ORR_event_rate,
    DCR_event_rate = DCR_event_rate,
    ORR=ORR,
    DCR=DCR
  )

  cat("\nORR(95%CI):\n")
  print(paste0(ORR_event_rate,ORR))
  cat("\nDCR(95%CI):\n")
  print(paste0(DCR_event_rate,DCR))
  return(result)
}

####卡方 计算####
#' Used to quickly calculate chisq or fisher values
#'
#' @title chisq_fisher
#' @description Performs Chi-squared or Fisher's exact test based on user-specified contingency table data.
#' Automatically detects whether to use Yates correction or Fisher's test.
#' @param data A numeric vector representing contingency table counts, filled row-wise.
#' @param ncol Integer. Number of columns in the contingency table.
#' @param row_names Optional character vector for row labels.
#' @param col_names Optional character vector for column labels.
#' @param language Output language. Either "en" (default) or "cn".
#' @author camellia
#' @examples
#' chisq_fisher(
#'   data = c(15, 25, 10, 13),
#'   ncol = 2
#' )
#' @return A list including the contingency table, expected values, p-value, method used and other metadata.
#' @export
chisq_fisher <- function(data,
                         ncol,
                         row_names = NULL,
                         col_names = NULL,
                         digits = 2,
                         p_digits = 3) {
  if (missing(data) || missing(ncol)) {
    stop("Both 'data' and 'ncol' are required.")
  }
  if (length(data) %% ncol != 0) {
    stop("Data length must be a multiple of 'ncol'.")
  }

  nrow <- length(data) / ncol
  mat <- matrix(data, nrow = nrow, ncol = ncol, byrow = TRUE)
  if (is.null(row_names)) row_names <- paste0("Row", 1:nrow)
  if (is.null(col_names)) col_names <- paste0("Col", 1:ncol)
  rownames(mat) <- row_names
  colnames(mat) <- col_names

  cat("\n--- Contingency Table ---\n")
  print(mat)

  test_chi <- suppressWarnings(chisq.test(mat, correct = FALSE))
  expected <- round(test_chi$expected, digits)
  cat("\n--- Expected Frequencies ---\n")
  print(expected)

  min_expected <- min(expected)
  total_n <- sum(mat)
  prop_low <- mean(expected < 5)
  is_2x2 <- all(dim(mat) == c(2, 2))

  if (is_2x2) {
    if (min_expected >= 5 && total_n >= 40) {
      method_msg <- "Using Pearson's Chi-squared Test..."
      result <- suppressWarnings(chisq.test(mat, correct = FALSE))
    } else if (min_expected >= 1 && min_expected < 5 && total_n >= 40) {
      method_msg <- "Using Yates-corrected Chi-squared Test (2x2)..."
      result <- suppressWarnings(chisq.test(mat, correct = TRUE))
    } else {
      method_msg <- "Using Fisher's Exact Test (low expected counts or small N)..."
      result <- suppressWarnings(fisher.test(mat))
    }
  } else {
    if (min_expected >= 1 && prop_low <= 0.2) {
      method_msg <- "Using Pearson's Chi-squared Test..."
      result <- suppressWarnings(chisq.test(mat, correct = FALSE))
    } else {
      method_msg <- "Using Fisher's Exact Test (low expected counts)..."
      result <- suppressWarnings(fisher.test(mat))
    }
  }

  stat_val <- if (is.numeric(result$statistic)) sprintf(paste0("%.0", digits, "f"), result$statistic) else result$statistic
  df_val <- if (!is.null(result$parameter)) sprintf(paste0("%.0", digits, "f"), result$parameter) else NA

  format_p <- function(p) {
    if (is.na(p)) return("p = NA")
    if (p < 10^(-p_digits)) {
      paste0("p < ", sprintf(paste0("%.0", p_digits, "f"), 10^(-p_digits)))
    } else {
      paste0("p = ", sprintf(paste0("%.0", p_digits, "f"), p))
    }
  }

  p_label <- format_p(result$p.value)

  cat("\n", method_msg, "\n")
  cat("\n--- Test Summary ---\n")
  print(result)
  cat("\n---p-value---\n")
  print(paste0(result$p.value,"    ",p_label))

  invisible(list(
    table = mat,
    expected = expected,
    method_used = result$method,
    method_description = method_msg,
    statistic = stat_val,
    df = df_val,
    p.value = result$p.value,
    p.label = p_label,
    result = result
  ))
}

#' Draw Group Comparison Plot with Statistical Tests
#'
#' This function generates a group comparison plot for multiple timepoints or variables.
#' It performs normality and variance homogeneity tests (Shapiro and Levene), then
#' applies appropriate statistical comparisons (t-test or Wilcoxon test), and summarizes results in a GT table.
#'
#' @param data A data frame containing the input data.
#' @param value_cols A character vector of column names representing different timepoints or measurements.
#' @param eff_col A character string indicating the column name for grouping (e.g., treatment group).
#' @param type_labels Optional named vector to relabel `value_cols` in the plot.
#' @param y_label Y-axis label for the plot. Default is `"B cell(%)"`.
#' @param y_limit A numeric vector of length 2 setting Y-axis limits. Default is `c(0, 4)`.
#' @param y_breaks A numeric vector specifying the breaks for Y-axis. Default is `seq(0, 4, 1)`.
#' @param palette A vector of colors used for plotting different groups. Default is a blue-yellow palette.
#' @param p_label_y Optional numeric value or vector to control vertical position of p-value labels.
#' @param dodge_width Width for dodging points and error bars. Default is 0.5.
#' @param p_digits Number of decimal places to show in p-values. Default is 3.
#' @return A list with the following components:
#' \item{plot}{A `ggplot2` object showing the comparison plot.}
#' \item{results}{A data frame of raw statistical test results.}
#' \item{formatted_results}{A formatted table of results with significance, normality, and variance assumptions.}
#' \item{gt_table}{A `gt` object showing a styled result table (if `gt` is installed).}
#'
#' @import dplyr ggplot2 tidyr broom car gt
#' @example
#' set.seed(123)
#'library(ggplot2)
#'library(gt)
#'library(tidyverse)
#'library(readxl)
#'library(broom)
#'library(car)
#' n <- 100
#' sim_data <- data.frame(
#'  受试者编号 = paste0("S", 1:n),
  #'  疗效 = rep(c("ORR", "Non-ORR"), each = n/2),
  #'  免疫诱导前 = c(
    #'   rnorm(n/2, mean = 60, sd = 10),
    #'   rnorm(n/2, mean = 50, sd = 15)
    #'  ),
  #'  免疫诱导后 = c(
    #'   rnorm(n/2, mean = 65, sd = 8),
    #'   rnorm(n/2, mean = 55, sd = 12)
    #'   )
  #' )
#' head(sim_data)
#' type_map <- c(
#'  "免疫诱导前" = "Pre",
#'  "免疫诱导后" = "Post"
#' )
#' result <- draw_group_comparison_plot(
#'  data = sim_data,
#'  value_cols = c(3:6),#指点列转长格式
#'  eff_col = "疗效",#分组变量
#'   type_labels = type_map,#x轴标签
#'  y_label = "T lymphocyte content",
#'  y_limit = c(10, 100),
#'  y_breaks = seq(10, 100, 20),
#'  p_label_y=c(90),#p值显示位置
#'  dodge_width=0.1,#偏倚
#'   palette = c("#1B9E77", "#D95F02"),
#'   p_digits = 3
#' )
#' result$plot
#' result$gt_table
#' @export
#---组间比较 绘图 输出p值----
draw_group_comparison_plot <- function(data, value_cols, eff_col,
                                       type_labels = NULL, y_label = "B cell(%)",
                                       y_limit = c(0, 4), y_breaks = seq(0, 4, 1),
                                       palette = c("#0073C2FF", "#EFC000FF"),
                                       p_label_y = NULL, dodge_width = 0.5,
                                       p_digits = 3) {
  dt_long <- data %>%
    pivot_longer(cols = all_of(value_cols),
                 names_to = "type", values_to = "value") %>%
    na.omit()

  if (!is.null(type_labels)) {
    dt_long$type <- factor(dt_long$type,
                           levels = names(type_labels), labels = type_labels)
  }

  norm_levene <- dt_long %>%
    group_by(type) %>%
    summarise(
      p_shapiro = shapiro.test(value)$p.value,
      p_levene = leveneTest(value ~ as.factor(get(eff_col)))[[1, "Pr(>F)"]],
      .groups = "drop"
    )

  data_split <- dt_long %>% group_split(type)

  test_fun <- function(df) {
    type_name <- unique(df$type)
    p_shapiro <- norm_levene$p_shapiro[norm_levene$type == type_name]
    p_levene <- norm_levene$p_levene[norm_levene$type == type_name]

    if (p_shapiro > 0.05) {
      if (p_levene > 0.05) {
        res <- t.test(value ~ get(eff_col), data = df, var.equal = TRUE)
        method <- "t.test (equal var)"
      } else {
        res <- t.test(value ~ get(eff_col), data = df)
        method <- "t.test (Welch)"
      }
    } else {
      res <- wilcox.test(value ~ get(eff_col), data = df)
      method <- "wilcox.test"
    }

    broom::tidy(res) %>%
      mutate(type = type_name,
             method = method,
             p_shapiro = p_shapiro,
             p_levene = p_levene)
  }

  results <- lapply(data_split, test_fun) %>%
    bind_rows() %>%
    mutate(
      label = ifelse(p.value < 0.001, "p < 0.001",
                     paste0("p = ", sprintf(paste0("%.", p_digits, "f"), p.value)))
    )

  if (is.null(p_label_y)) {
    y_range <- diff(range(dt_long$value, na.rm = TRUE))
    results$y_pos <- max(dt_long$value, na.rm = TRUE) + 0.1 * y_range
  } else {
    if (length(p_label_y) == 1) {
      results$y_pos <- p_label_y
    } else if (length(p_label_y) == nrow(results)) {
      results$y_pos <- p_label_y
    } else {
      warning("p_label_y长度不匹配，使用默认位置")
      y_range <- diff(range(dt_long$value, na.rm = TRUE))
      results$y_pos <- max(dt_long$value, na.rm = TRUE) + 0.1 * y_range
    }
  }

  format_p_value <- function(p) {
    ifelse(p < 0.001, "<0.001", sprintf(paste0("%.", p_digits, "f"), p))
  }

  formatted_results <- results %>%
    select(type, method, p_shapiro, p_levene, p.value, estimate) %>%
    mutate(
      p_shapiro = format_p_value(p_shapiro),
      p_levene = format_p_value(p_levene),
      p.value = format_p_value(p.value),
      estimate = sprintf("%.4f", estimate),
      `正态性检验` = ifelse(p_shapiro > 0.05, "正态", "非正态"),
      `方差齐性` = ifelse(p_levene > 0.05, "齐性", "不齐"),
      `显著性` = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        TRUE ~ "不显著"
      )
    ) %>%
    select(
      `时间点` = type,
      `检验方法` = method,
      `正态性检验p值` = p_shapiro,
      `正态性` = `正态性检验`,
      `方差齐性p值` = p_levene,
      `方差齐性` = `方差齐性`,
      `组间比较p值` = p.value,
      `显著性` = `显著性`,
      `效应量估计` = estimate
    )

  if (requireNamespace("gt", quietly = TRUE)) {
    gt_table <- formatted_results %>%
      gt::gt() %>%
      gt::tab_header(
        title = "组间比较统计结果",
        subtitle = paste0(eff_col, "组在不同时间点的比较")
      ) %>%
      gt::tab_style(
        style = gt::cell_text(weight = "bold"),
        locations = gt::cells_column_labels()
      ) %>%
      gt::tab_style(
        style = gt::cell_fill(color = "lightgreen"),
        locations = gt::cells_body(
          columns = "显著性",
          rows = `显著性` == "***"
        )
      ) %>%
      gt::tab_style(
        style = gt::cell_fill(color = "lightyellow"),
        locations = gt::cells_body(
          columns = "显著性",
          rows = `显著性` == "**"
        )
      ) %>%
      gt::tab_style(
        style = gt::cell_fill(color = "lightpink"),
        locations = gt::cells_body(
          columns = "显著性",
          rows = `显著性` == "*"
        )
      )
  } else {
    warning("gt包未安装，无法创建美化表格")
    gt_table <- NULL
  }

  summary_data <- dt_long %>%
    group_by(type, !!sym(eff_col)) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      se = sd(value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )

  p <- ggplot(dt_long, aes(x = type, y = value, color = .data[[eff_col]])) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2,
                                                dodge.width = dodge_width),
                alpha = 0.6, size = 1) +
    geom_errorbar(data = summary_data,
                  aes(y = mean, ymin = mean - se, ymax = mean + se,
                      color = .data[[eff_col]]),
                  width = 0.1, size = 0.8,
                  position = position_dodge(width = dodge_width)) +
    geom_point(data = summary_data,
               aes(y = mean, fill = .data[[eff_col]]),
               size = 2.5, shape = 21, color = "black",
               position = position_dodge(width = dodge_width)) +
    geom_text(data = results,
              aes(x = type, y = y_pos, label = label),
              inherit.aes = FALSE,
              size = 4) +
    geom_line(data = summary_data,
              aes(x = type, y = mean, group = !!sym(eff_col),
                  color = !!sym(eff_col)),
              size = 0.8, position = position_dodge(width = dodge_width)) +
    scale_y_continuous(breaks = y_breaks, limits = y_limit, expand = c(0.1, 0.1)) +
    scale_color_manual(values = palette) +
    scale_fill_manual(values = palette) +
    theme_classic(base_size = 14) +
    labs(y = y_label, x = NULL, color = NULL, fill = NULL,
         caption = "Error bars represent standard error of the mean (SEM).") +
    theme(legend.position = "top",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          plot.caption = element_text(size = 10))

  return(list(
    plot = p,
    results = results,
    formatted_results = formatted_results,
    gt_table = gt_table
  ))
}
