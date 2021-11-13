library(dplyr)
library(sandwich)
library(lmtest)

# для проверки мощности ставить group > 0
true_coefs <- list(
  'intercept' = 7, 
  'group' = 0.3, 
  'network' = 2
  )

results <- list(
  'treated_share' = c(),
  'group_estimated' = c()
  )

boxplotData <- data.frame()

#читаем директорию со сгенерированными графами
files  <- list.files("~/linear-models-on-graphs/generated_graphs/")

for (graph_file in files) {
  df <- read.csv2(
    paste0('~/linear-models-on-graphs/generated_graphs/', graph_file), 
    sep = ',', 
    dec = '.'
  )
  
  # считаем показатели по всем вершинам
  treated_friends_cnt <- df %>% 
    group_by(vertex_1) %>% 
    summarise(
      treated_share = mean(is_treated),
      friends = length(vertex_2),
      treated_friends = sum(is_treated==1),
      group_id = ifelse(unique(vertex_1) > 49, 1, 0)
    )
  
  # доля тритед друзей в тестовых группах
  treated_shares <- treated_friends_cnt %>% 
    group_by(group_id) %>%
    summarise(treated_avg = mean(treated_share))
  
  # генерируем таргет
  treated_friends_cnt <- treated_friends_cnt %>%
    mutate(
      target = true_coefs$intercept + true_coefs$group * group_id + true_coefs$network * treated_friends + rnorm(100, 0, 1),
      cluster_id = treated_friends # sample(0:1, 100, replace = TRUE) #ifelse(unique(vertex_1) > 49, 1, 0)
    )
  
  result_lm <- lm(
    target ~ group_id + treated_friends, 
    data = treated_friends_cnt
  )
  estimated_coefs <- result_lm$coefficients
  model_summary <- summary(result_lm)
  
  # вспомогательная регрессия для подсчета RSS-j
  additional_lm <- lm(
    group_id ~ treated_friends, 
    data = treated_friends_cnt
  )

  # значение коэффициента при группе, величина модулярности и доля тритед друзей в контроле
  key <- substr(graph_file, 18, 20)
  boxplotData <- bind_rows(
    boxplotData, 
    data.frame(
      modularity = key, 
      coef = estimated_coefs[2],
      stderr = model_summary$coefficients[2, c('Std. Error')],
      stderr_nw = sqrt(NeweyWest(result_lm)[2,2]),
      stderr_cl = coeftest(result_lm, vcov = vcovCL, cluster = ~cluster_id)[2,2],
      rss = sum(result_lm$residuals^2)/(nrow(treated_friends_cnt)-length(true_coefs)-1),
      rss_j = sum(additional_lm$residuals^2)/(nrow(treated_friends_cnt)-length(true_coefs)-2),
      var_treated_grp = var(treated_friends_cnt$treated_friends[1:50])/var(treated_friends_cnt$treated_friends[51:100]),
      var_treated = var(treated_friends_cnt$treated_friends),
      bpt_pval= lmtest::bptest(result_lm)$p.value,
      dwt_pval = lmtest::dwtest(result_lm)$p.value,
      swt_pval = shapiro.test(result_lm$residuals)$p.value,
      treated_share = treated_shares$treated_avg[1],
      pval_lm = model_summary$coefficients[2,c('Pr(>|t|)')],
      pval_nw = 2* pt(
        abs(estimated_coefs[2]/sqrt(NeweyWest(result_lm)[2,2])), 
        result_lm$df.residual, 
        lower.tail = F
        ),
      pval_cl = 2* pt(
        abs(estimated_coefs[2]/coeftest(result_lm, vcov = vcovCL, cluster = ~cluster_id)[2,2]), 
        result_lm$df.residual, 
        lower.tail = F
        )
      )
    )
}

# смотрим стандартные ошибки 
boxplotData %>%
  group_by(modularity) %>%
  summarise(
    avg_stderr = mean(stderr),
    avg_stderr_newey_west = mean(stderr_nw),
    avg_stderr_clustered = mean(stderr_cl),
  )

# смотрим диагностику справедливости классических предпосылок
boxplotData %>%
  group_by(modularity) %>%
  summarise(
    BreuschPaganTest = sum(bpt_pval < 0.05)/n() * 100,
    DurbinWatsonTest = sum(dwt_pval < 0.05)/n() * 100,
    ShapiroWilkTest = sum(swt_pval < 0.05)/n() * 100
  )

# ошибки первого и второго рода
boxplotData %>%
  group_by(modularity) %>%
  summarise(
    pow_classic = sum(pval_lm < 0.05)/n() * 100,
    pow_newey_west = sum(pval_nw < 0.05)/n() * 100,
    pow_clustered = sum(pval_cl < 0.05)/n() * 100
  )

### порисуем графики 

plot(
  boxplotData$treated_share, 
  boxplotData$coef, 
  main="treated share vs coef estimation",
  xlab="treated_share", 
  ylab="group_estimated", 
  pch=1
  )

boxplot(
  coef ~ modularity, 
  data = boxplotData,
  main="generated modularity and group estimation"
  )

boxplot(
  stderr ~ modularity, 
  data = boxplotData,
  main="generated modularity and group stderr"
)

boxplot(
  stderr_nw ~ modularity, 
  data = boxplotData,
  main="generated modularity and group stderr-NW"
)

boxplot(
  stderr_cl ~ modularity, 
  data = boxplotData,
  main="generated modularity and group stderr-clustered"
)

boxplot(
  rss ~ modularity, 
  data = boxplotData,
  main="generated modularity and RSS"
)

boxplot(
  rss_j ~ modularity, 
  data = boxplotData,
  main="generated modularity and RSS-j"
)

boxplot(
  treated_share ~ modularity, 
  data = boxplotData,
  main="generated modularity and treated_share"
)

boxplot(
  var_treated ~ modularity, 
  data = boxplotData,
  main="generated modularity and var_treated"
)

plot(
  treated_friends_cnt$treated_friends, 
  treated_friends_cnt$target, 
  main="treated_friends vs target",
  xlab="treated_friends", 
  ylab="target", 
  pch=1
)

# looking at residuals
plot(
  treated_friends_cnt$target,
  result_lm$residuals
)

plot(result_lm$residuals[order(result_lm$residuals)])

