

x <- c(1:10)
class(x)<- "sexy"

print.sexy <-function(variables) {
  print("im really sexy")

}
summary.sexy <- function(Variable){
  print("to sexy for all these cloths")
}
ts.sexy <- function(Variable){
  print("to sexy for sleep")
}
print(x)
summary(x)
ts(x)
