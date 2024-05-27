G_2015 <- daymeann2015
G_2016 <- daymeann2016
G_2017 <- daymeann2017
G_2018 <- daymeann2018
G_2019 <- daymeann2019
G_2020 <- daymeann2020

K_2015 <- daymean2015
K_2016 <- daymean2016
K_2017 <- daymean2017
K_2018 <- daymean2018
K_2019 <- daymean2019
K_2020 <- daymean2020

#hypothesis: water level is the same
# if p-value is small, there is a significant difference between the two datasets 

#2015
# p-value: 2.2e-16 
t.test(G_2015$Average, K_2015$Average)

#2016
# p-value: 2.2e-16 
t.test(G_2016$Average, K_2016$Average,var.eq = T)

#2017 - too little data
# p-value: 2.2e-16 
t.test(G_2017$Average, K_2017$Average,var.eq = T)

#2018
# p-value: 2.2e-16 
t.test(G_2018$Average, K_2018$Average,var.eq = T)

#2019
# p-value: 0.05627 
t.test(G_2019$Average, K_2019$Average,var.eq = T)

#2020 - too little data
# p-value: 2.2e-16 
t.test(G_2020$Average, K_2020$Average,var.eq = T)
