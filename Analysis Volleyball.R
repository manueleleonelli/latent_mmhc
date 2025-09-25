setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Current Papers/BNs Volleyball")

volley_data <- as.data.frame(readRDS("volley_data.rds"))

bigfive <- c("neuroticism_avg","agreeableness_avg","conscientiousness_avg","extraversion_avg","openness_avg")
psycho <- c("self_talk_avg", "goal_setting_avg", "self_conf_avg", "emotional_arousal_avg", "worry_avg", "concentration_disruption_avg", "mental_practice_avg", "match_preparation_avg")
self_esteem <- c("self_esteem_avg")
var <- c("ruolgioc","serie","altsport","oreallapv","compuff","arbitrsport","condatt","relaz","fumo","bevalcol","Eta")

data <- volley_data[,c(bigfive,psycho, self_esteem, var)]

data$ruolgioc <- factor(data$ruolgioc)
data$serie  <- factor(data$serie)
data$altsport <- factor(data$altsport)
data$oreallapv <- factor(data$oreallapv,ordered=T)
data$compuff <- factor(data$compuff, ordered = T)
data$arbitrsport <- factor(data$arbitrsport)
data$condatt <- factor(data$condatt)
data$relaz <- factor(data$relaz)
data$fumo <- factor(data$fumo)
data$bevalcol <- factor(data$bevalcol)
data$Eta <- cut(data$Eta,breaks = c(18, 20, 25, 30, Inf), labels = c("18–19", "20–24", "25–29", "30+"), right = FALSE, ordered_result = TRUE)

data$ruolgioc <- factor( ifelse(data$ruolgioc %in% c("Opposite Hitter", "Wing Spiker", "Middle Blocker"),"Hitter", "Non_Hitter"),levels = c("Non_Hitter", "Hitter"))

data$condatt <- factor(ifelse(data$condatt %in% c("Worker", "Both"),"Working", "Not_Working"),levels = c("Not_Working", "Working"))



# Indices
non_demo_idx <- 1:14
demo_idx     <- 15:25

# Names
vn <- colnames(data)
from_vars <- vn[non_demo_idx]
to_vars   <- vn[demo_idx]

# Blacklist: forbid edges FROM 1:14 TO 15:25
blacklist_demo <- expand.grid(
  from = from_vars,
  to   = to_vars,
  stringsAsFactors = FALSE
)
V  <- colnames(data)
p  <- length(V)
n  <- nrow(data)

fixedGaps <- matrix(FALSE, p, p, dimnames = list(V, V))

non_demo <- V[1:14]
demo     <- V[15:25]

# forbid any edge between {1..14} and {15..25}
fixedGaps[non_demo, demo] <- TRUE
fixedGaps[demo, non_demo] <- TRUE  # symmetric!


# Example: use in bnlearn learners (hc/tabu/mmhc)
# dag <- bnlearn::tabu(data, score = "bic", blacklist = blacklist_demo)

# If you already have a blacklist (e.g., from a PC skeleton step), union them:
# blacklist <- unique(rbind(blacklist, blacklist_demo))

learn_dag_mmhc_sem <- function(data, corr, alpha = 0.01,
                               use_cholesky = TRUE,
                               verbose = FALSE,
                               penalty_type = "BIC",
                               blacklist = NULL) {
  require(pcalg)
  require(bnlearn)
  
  stopifnot(is.data.frame(data))
  var.names <- colnames(data)
  n <- nrow(data)
  
  # PC skeleton from corr
  suffStat <- list(C = corr, n = n)
  pc_fit <- pc(suffStat = suffStat, indepTest = gaussCItest,
               labels = var.names, alpha = alpha, conservative = TRUE, fixedGaps = fixedGaps)
  amat_pc <- as(pc_fit@graph, "matrix")
  
  # skeleton -> blacklist
  skeleton <- which((amat_pc + t(amat_pc)) > 0, arr.ind = TRUE)
  skeleton <- skeleton[skeleton[,1] != skeleton[,2], , drop = FALSE]
  skel_names <- data.frame(from = var.names[skeleton[,1]],
                           to   = var.names[skeleton[,2]], 
                           stringsAsFactors = FALSE)
  all_edges <- expand.grid(from = var.names, to = var.names, stringsAsFactors = FALSE)
  all_edges <- all_edges[all_edges$from != all_edges$to, ]
  keep_edges <- rbind(skel_names, skel_names[, c("to","from")])
  edge_str <- paste(keep_edges$from, keep_edges$to, sep = "->")
  all_str  <- paste(all_edges$from, all_edges$to,  sep = "->")
  blacklist <- all_edges[!all_str %in% edge_str, ]
  
  # latent Gaussian Z from ranks; optional Cholesky whiten
  ranks <- apply(data, 2, rank, ties.method = "average") / (n + 1)
  Z_raw <- qnorm(ranks)
  if (use_cholesky) {
    chol_corr <- tryCatch(chol(corr), error = function(e) {
      eps <- 1e-6; chol((1 - eps) * corr + eps * diag(ncol(corr)))
    })
    Z <- scale(Z_raw, center = TRUE, scale = FALSE) %*% chol_corr
  } else {
    Z <- Z_raw
  }
  
  # pass penalty_type here
  args_list <- list(Z = Z, n = n, penalty_type = penalty_type)
  
  bnlearn::tabu(
    data,
    score = "custom-score",
    fun = sem_copula_score,   # <-- 4-arg function above
    args = args_list,
    blacklist = blacklist
  )
}

cop_ken <- estimate_copula_correlation_kendall(data)
mmhc_sem_05_raw = learn_dag_mmhc_sem(data, cop_ken, alpha = 0.01, use_cholesky = FALSE,blacklist=blacklist_demo)
bot <- learn_dag_mmhc_sem_bootstrap(data, B = 500, alpha = 0.10,blacklist=blacklist_demo)


data <- volley_data[,c(bigfive,psycho, self_esteem)]
data$neuroticism <- factor(round(data$neuroticism_avg),levels = 1:5,ordered = T)
data$agreeableness<- factor(round(data$agreeableness_avg),levels = 1:5,ordered = T)
data$conscientiousness <- factor(round(data$conscientiousness_avg),levels = 1:5,ordered = T)
data$extraversion <- factor(round(data$extraversion_avg),levels = 1:5,ordered = T)
data$openness <- factor(round(data$openness_avg),levels = 1:5,ordered = T)
data$self_talk <- factor(round(data$self_talk_avg),ordered=T)
data$goal_setting <- factor(round(data$goal_setting_avg), ordered = T)
data$self_conf <- factor(round(data$self_conf_avg), ordered = T)
data$emotional_arousal <- factor(round(data$emotional_arousal_avg), ordered = T)
data$worry <- factor(round(data$worry_avg), ordered = T)
data$concentration_disruption <- factor(round(data$concentration_disruption_avg), ordered = T)
data$concentration_disruption <- factor(data$concentration_disruption,levels = 1:6,ordered = TRUE)
data$mental_practice <- factor(round(data$mental_practice_avg), ordered = T)
data$match_preparation <- factor(round(data$match_preparation_avg), ordered = T)
data$self_esteem <- factor(round(data$self_esteem_avg), ordered = T)
data$self_esteem <- factor(data$self_esteem,levels = 1:4,ordered = TRUE)
data <- data[,-c(1:14)]

# Original arcs as a matrix
arcs <- matrix(c(
  "neuroticism",              "conscientiousness",
  "neuroticism",              "worry",
  "conscientiousness",        "agreeableness",
  "conscientiousness",        "openness",
  "extraversion",             "agreeableness",
  "extraversion",             "self_esteem",
  "goal_setting",             "self_conf",
  "goal_setting",             "mental_practice",
  "goal_setting",             "match_preparation",
  "self_conf",                "mental_practice",
  "self_conf",                "self_esteem",
  "emotional_arousal",        "extraversion",
  "emotional_arousal",        "goal_setting",
  "emotional_arousal",        "self_conf",
  "emotional_arousal",        "worry",
  "concentration_disruption", "worry",
  "mental_practice",          "match_preparation",
  "match_preparation",        "self_talk"
), byrow = TRUE, ncol = 2)

colnames(arcs) <- c("from", "to")

dag <- empty.graph(colnames(data))
arcs(dag) <- arcs
bn <- bn.fit(dag,data,method="bayes", iss=1)

write.net("volley.net",bn)

names(data) <- c(
  "player_position", "league_division", "plays_other_sports", "weekly_training_hours",
  "num_competitions", "referee_experience", "working_status", "in_relationship",
  "smokes", "drinks_alcohol", "age_group",
  "neuroticism", "agreeableness", "conscientiousness", "extraversion", "openness",
  "self_talk", "goal_setting", "self_confidence", "emotional_arousal", "worry",
  "concentration_issues", "mental_practice", "match_preparation", "self_esteem"
)
