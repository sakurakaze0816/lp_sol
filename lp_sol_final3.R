# 誤差容忍值
min_weight_default <- 0.000001

# 辨認與安裝
toInstall <- c("lpSolve","openxlsx","sqldf","magrittr")
packages <- installed.packages()
whichInstall <- which(!(toInstall %in% packages))
for(i in whichInstall){
  istl <- toInstall[i]
  print(paste("*** 系統安裝套件 (第一次執行時)：",istl,sep=""))
  install.packages(istl, quiet = TRUE , repos="https://cloud.r-project.org")
  print(paste("*** 成功安裝：",istl,sep=""))
}
print("*** 開始：載入套件")

library(sqldf)
library(lpSolve)
library(openxlsx)

#讀檔
lp_data <- read.xlsx(xlsxFile="lp_data.xlsx", sheet = 1, skipEmptyRows = FALSE, colNames = TRUE, na.strings = "")
cc <- cbind(lp_data$mean_norm,lp_data$TSignal_norm, lp_data$wc_norm, lp_data$wp_norm)

lp_data <- lp_data[!duplicated(cc), ]                                                 # 去除重複資料
cc <- cbind(lp_data$mean_norm,lp_data$TSignal_norm, lp_data$wc_norm, lp_data$wp_norm) # alternatives

# start LP settings
f.con0 <- matrix(c(rbind(cc[,1],cc[,2],cc[,3],cc[,4])), nrow=nrow(cc), byrow=TRUE)    # eq(2)
f.dir0 <- c(rep("<=",nrow(cc)))
f.rhs0 <- rep(1,nrow(cc))

f.con <- matrix(c(rbind(cc[,1],cc[,2],cc[,3],cc[,4]),1,1,1,1), nrow=nrow(cc)+1, byrow=TRUE)   # eq(2) & eq(6)
f.dir <- c(rep("<=",nrow(cc)),'=')
f.rhs <- rep(1,nrow(cc)+1)

f.con1 <- f.con
f.dir1 <- f.dir
f.rhs1 <- f.rhs

for(i in 1:4){                                   # eq(3)
  row_e <- rep(0,4)
  row_e[i] <- 1
  f.con1 <- rbind(f.con1,row_e)
  f.dir1 <- c(f.dir1,">=")
  f.rhs1 <- c(f.rhs1,min_weight_default)
}

f.con2 <- f.con1
f.dir2 <- f.dir1
f.rhs2 <- f.rhs1

for(i in 1:3){                                  # eq(5)
  row_e <- rep(0,4)
  row_e[i] <- 1
  row_e[i+1] <- -1
  f.con2 <- rbind(f.con2,row_e)
  f.dir2 <- c(f.dir2,">=")
  f.rhs2 <- c(f.rhs2,0)
}

# 設定變數
objvals0 <- rep(0.0,nrow(cc))
objvals1 <- rep(0.0,nrow(cc))
objvals2 <- rep(0.0,nrow(cc))

solutions0 <- matrix(rep(0.0,4*nrow(cc)),nrow=nrow(cc),byrow=TRUE)
solutions1 <- matrix(rep(0.0,4*nrow(cc)),nrow=nrow(cc),byrow=TRUE)
solutions2 <- matrix(rep(0.0,4*nrow(cc)),nrow=nrow(cc),byrow=TRUE)

sijp1 <- matrix(rep(0.0,nrow(cc)*nrow(cc)),ncol=nrow(cc),nrow=nrow(cc),byrow=TRUE)
sijp2 <- matrix(rep(0.0,nrow(cc)*nrow(cc)),ncol=nrow(cc),nrow=nrow(cc),byrow=TRUE)

sip1 <- rep(0.0,nrow(cc))
sip2 <- rep(0.0,nrow(cc))

feasible0 <- rep('',nrow(cc))
feasible1 <- rep('',nrow(cc))
feasible2 <- rep('',nrow(cc))

pb   <- txtProgressBar(max = nrow(cc), style=3)

# solve LP

# eq(2)
# 求rank0的解
start_model_time <- Sys.time()                   # 程式運行時間
for(i in 1:nrow(cc)){                            # eq(4)
  f.obj0 <- c(cc[i,1],cc[i,2],cc[i,3],cc[i,4])   # for each item
  result <- lp("max",f.obj0,f.con0,f.dir0,f.rhs0)
  objvals0[i] <- result$objval
  solutions0[i,] <- result$solution
  feasible0[i] <- if(result$status == 2){        # 判斷是否可行
    "不可行"
  }else{
    "可行"
  }
  setTxtProgressBar(pb, i)                       # 程式運行進度條
}
end_model_time <- Sys.time()
duration0 <- as.numeric(difftime(end_model_time, start_model_time, units = "secs"))

# eq(2) & eq(3) & eq(6)
# 求rank1的解
start_model_time <- Sys.time()                   # 程式運行時間
for(i in 1:nrow(cc)){                            # eq(4)
  f.obj1 <- c(cc[i,1],cc[i,2],cc[i,3],cc[i,4])   # for each item
  result <- lp("max",f.obj1,f.con1,f.dir1,f.rhs1)
  objvals1[i] <- result$objval
  solutions1[i,] <- result$solution
  feasible1[i] <- if(result$status == 2){        # 判斷是否可行
    "不可行"
  }else{
    "可行"
  }
  setTxtProgressBar(pb, i)                       # 程式運行進度條
}
end_model_time <- Sys.time()
duration1 <- as.numeric(difftime(end_model_time, start_model_time, units = "secs"))

# eq(2) & eq(3) & eq(6)
# 求rank2的解
start_model_time <- Sys.time()                   # 程式運行時間
for(i in 1:nrow(cc)){
  for(j in 1:nrow(cc)){
    sijp1[j,i] <- sum(solutions1[j,]*cc[i,])     # eq(4')
  }
  sip1[i] <- mean(sijp1[,i])                     # eq(5')
  setTxtProgressBar(pb, i)                       # 程式運行進度條
}
end_model_time <- Sys.time()
duration2 <- as.numeric(difftime(end_model_time, start_model_time, units = "secs"))

# eq(2) & eq(3) & eq(5) & eq(6)
# 求rank3的解
start_model_time <- Sys.time()                   # 程式運行時間
for(i in 1:nrow(cc)){                            # eq(4)
  f.obj2 <- c(cc[i,1],cc[i,2],cc[i,3],cc[i,4])   # for each item
  result <- lp("max",f.obj2,f.con2,f.dir2,f.rhs2)
  objvals2[i] <- result$objval
  solutions2[i,] <- result$solution
  feasible2[i] <- if(result$status == 2){        # 判斷是否可行
    "不可行"
  }else{
    "可行"
  }
  setTxtProgressBar(pb, i)                       # 程式運行進度條
}
end_model_time <- Sys.time()
duration3 <- as.numeric(difftime(end_model_time, start_model_time, units = "secs"))

# eq(2) & eq(3) & eq(5) & eq(6)
# 求rank4的解
start_model_time <- Sys.time()                   # 程式運行時間
for(i in 1:nrow(cc)){                                             
  for(j in 1:nrow(cc)){
    sijp2[j,i] <- sum(solutions2[j,]*cc[i,])     # eq(4')
  }
  sip2[i] <- mean(sijp2[,i])                     # eq(5')
  setTxtProgressBar(pb, i)                       # 程式運行進度條
}
end_model_time <- Sys.time()
duration4 <- as.numeric(difftime(end_model_time, start_model_time, units = "secs"))

#rank的變數設定
rank0 <- rank(-objvals0,ties.method="min")
rank1 <- rank(-objvals1,ties.method="min")
rank2 <- rank(-sip1,ties.method="min") 
rank3 <- rank(-objvals2,ties.method="min")
rank4 <- rank(-sip2,ties.method="min")

# 重複的數量
duplicate0 <- sum(duplicated(rank0))
duplicate1 <- sum(duplicated(rank1))        
duplicate2 <- sum(duplicated(rank2))
duplicate3 <- sum(duplicated(rank3))
duplicate4 <- sum(duplicated(rank4))

#合併所求的資料
lp_sol_final1 <- cbind(lp_data,objvals0,rank0,solutions0,feasible0,objvals1,rank1,solutions1,feasible1,sip1,rank2,objvals2,rank3,solutions2,feasible2,sip2,rank4)

#判斷重複的數量與執行的時間
cat(paste("rank0重複的數量",duplicate0,"執行的時間",duration0),"\n",
    paste("rank1重複的數量",duplicate1,"執行的時間",duration1),"\n",
    paste("rank2重複的數量",duplicate2,"執行的時間",duration2),"\n",
    paste("rank3重複的數量",duplicate3,"執行的時間",duration3),"\n",
    paste("rank4重複的數量",duplicate4,"執行的時間",duration4),sep="")

#寫入檔案
write.xlsx(lp_sol_final3,"lp_sol_final3.xlsx")
