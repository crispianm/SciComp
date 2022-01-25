data = [2, 3, 1, -1]

# for value in data: print(value**2)

def squared(list):
    squared_list = []

    for value in list:
        squared_list.append(value**2)

    return squared_list

def filter(list):
    squared_list = []

    for value in list:
        if value>=0:
            squared_list.append(value)

    return squared_list


# for value in data[::-1]: print(value)
    
def fact(x):
    numbers = range(x)
    factorial = 1
    for num in numbers:
        num += 1
        factorial *= num
    return factorial
    
print(fact(30))
import math
print(math.factorial(30))

# squared_list = filter(data)
# print(squared_list)