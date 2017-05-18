library(SuperExactTest)

testA = letters[1:10]
testB = letters[7:20]
testC = letters[10:25]

meta_list = list(testA, testB, testC)

names(meta_list) = c("testA", "testB", "testC")

meta_result = supertest(meta_list, n = 100) 

summary(meta_result)
