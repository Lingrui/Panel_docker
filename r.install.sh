R --vanilla <<EOF
#install.packages('Cairo',repos="http://mirrors.opencas.cn/cran/")
install.packages('Cairo', dependencies = TRUE, repos='http://cran.rstudio.com/')
q()
EOF
