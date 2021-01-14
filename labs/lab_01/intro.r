# Install the fda library & call the fda library
if (!require("pacman")) install.packages("pacman")
pacman::p_load(fda, psych)
# Description of Canadian Weather data
data(CanadianWeather)
summary(CanadianWeather)
# Extract the temperature data of the 35 cities
temp = CanadianWeather$dailyAv
df_1 = as.data.frame(CanadianWeather$dailyAv)[,1:35]
# And plot them using matplot 
matplot(df_1,type="l",xlab="Days")
# Extract the Vancouver temperature data and plot it
plot(df_1$Vancouver.Temperature.C,xlab="Days")
# Extract the Vancouver precipitation data and plot it
df_2 = as.data.frame(CanadianWeather$dailyAv)[,36:70]
plot(df_2$Vancouver.Precipitation.mm,xlab="Days")
# write out data to CSV
df = as.data.frame(CanadianWeather$dailyAv)
write.csv(df, "/cloud/project/labs/lab_01/CanadianWeather.csv")
