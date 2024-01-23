# This Python program demonstrates basic input, calculation, and output operations.

# Print a simple greeting to the console
print("Hello world")
# Read an integer input from the user and store it in the variable num1
num1 = int(input("Enter the first number: "))

# Read another integer input from the user and store it in the variable num2
num2 = int(input("Enter the second number: "))

# Calculate the sum of the two input numbers
sum1 = num1 + num2 

# Print the result along with an explanatory message
print("You have now recorded two numbers whose sum is:", sum1)

# Using an f-string to format the output message with the calculated sum
print(f"You have now recorded two numbers whose sum is: {sum1}")

# Creating a message by concatenating strings and then printing it
message = "You have now recorded two numbers whose sum is: " + str(sum1)
print(message)

# Using string formatting with the `.format()` method to insert the sum into the message
message = "You have now recorded two numbers whose sum is: {}".format(sum1)
print(message)
# Input two numbers from the user
num1 = int(input("Enter the first number: "))
num2 = int(input("Enter the second number: "))

# Calculate the sum of the two numbers
sum1 = num1 + num2
x = 10
y = False
# Check if the sum is greater than 10
x = int(input("print x"))
if sum1 > 10:
    print("The sum is greater than 10.")
    print(30)
    if x > 200:
        print(200)
elif sum1 > 20  and ( x > 100 or x > 90) :
    print(30)
    print(20)
else:
    print("The sum is not greater than 10.")


    # Input two numbers from the user
num1 = int(input("Enter the first number: "))
num2 = int(input("Enter the second number: "))

# Calculate the sum of the two numbers
sum1 = num1 + num2

# Check if the sum is greater than 10
if sum1 > 10:
    # If the sum is greater than 10, this block of code will execute
    print("The sum is greater than 10.")
else:
    # If the sum is not greater than 10, this block of code will execute
    print("The sum is not greater than 10.")
    
# The program execution continues here after the if-else block
for num in range(1, 6):
    print(num)

sum = 0
for num in range(1, 11):
    sum += num
print("Sum:", sum)

word = "Python" # [p, y, t, h, o , n]
for char in word:
    print(char)

for i in range(1, 6):
    for j in range(1, 11):
        print(i, "*", j, "=", i * j)

# The enumerate function can be used to iterate over both 
# the index and the value of a sequence.
fruits = ["apple", "banana", "cherry"]
for index, fruit in enumerate(fruits):
    print("Index:", index, "Fruit:", fruit)

for item in fruits:
    print(item)
num = 0
if num > 10:
    while num <= 1:
        print(num)
        num += 0.1
# Basic string creation
greeting = "Hello, World!"

# String concatenation
name = "Alice"
message = greeting + " My name is " + name
print(message)

# String indexing and slicing
print("First character:", greeting[0])
print("Last character:", greeting[-1])
print("Characters from index 7 to end:", greeting[7:])
print("Substring:", greeting[7:12])

# String length
length = len(greeting)# greeting.length 
print("Length of the string:", length)

# String methods
uppercase = greeting.upper()
lowercase = greeting.lower()
capitalized = greeting.capitalize()
replaced = greeting.replace("Hello", "Hi")

print("Uppercase:", uppercase)
print("Lowercase:", lowercase)
print("Capitalized:", capitalized)
print("Replaced:", replaced)

# String formatting
age = 25
formatted_string = "I am {} years old.".format(age)
print(formatted_string)

# f-string formatting (Python 3.6+)
formatted_string_f = f"I am {age} years old."
print(formatted_string_f)

# String splitting and joining
sentence = "This,is,a,sentence, to,split."
words = sentence.split(",")# [This,is,a,is.....   ]
joined_sentence = "-".join(words)
print("Words:", words)
print("Joined:", joined_sentence)
# In these examples:
# Slicing is done using the syntax [start:end], where start is the index of the first character you want to include, and end is the index after the last character you want to include.
# Omitting start means the slice starts from the beginning of the string.
# Omitting end means the slice goes until the end of the string.
# Negative indices count from the end of the string.
# A step value can be added as a third parameter to slice with a certain interval.
# Reversing a string can be achieved by using a negative step value.
# Remember that Python uses 0-based indexing, so the first character is at index 0, the second at index 1, and so on.

text = "Python Programming"

# Slicing from index 0 to index 5 (excluding index 5)
substring1 = text[0:5]
print(substring1)  # Output: "Python"

# Slicing from index 7 to the end of the string
substring2 = text[7:]
print(substring2)  # Output: "Programming"

# Slicing from index 7 to index 15 (excluding index 15)
substring3 = text[7:15]
print(substring3)  # Output: "Programm"

# Slicing with negative indices (counting from the end of the string)
substring4 = text[-11:-1]
print(substring4)  # Output: "Programmin"

# Slicing with a step (every second character)
substring5 = text[0::2]
print(substring5)  # Output: "Pto rgamn"

# Reversing a string using slicing
reversed_text = text[::-1]
print(reversed_text)  # Output: "gnimmargorP nohtyP"
#There are two common uses of triple quotes:

#1.Single-Line Strings: You can use triple quotes to create single-line strings just like regular quotes:
single_line_string = """This is a single-line string."""

# 2.Multi-Line Strings (Docstrings):
# Triple quotes are often used to create multi-line strings
# for documenting functions, classes, or modules.
#  This is called a docstring. It's a common practice to include 
# docstrings at the beginning of these constructs to explain their 
# purpose and usage:

def my_function(arg1, arg2):
    """
    This is a docstring that explains what the function does.
    
    Args:
        arg1 (int): The first argument.
        arg2 (str): The second argument.
        
    Returns:
        bool: The result of some operation.
    """
    # Function code here
def greet(name):
    """This function greets the person passed in as a parameter."""
    print(f"Hello, {name}! How are you?")

def add(a, b):
    """This function adds two numbers and returns the result."""
    result = a + b
    return result

def calculate_circle_area(radius):
    """This function calculates the area of a circle given its radius."""
    pi = 3.14159
    area = pi * radius ** 2
    return area


# Calling the greet function
greet("Alice")
greet("Bob")

# Calling the add function
sum_result = add(5, 3)
print("Sum:", sum_result)

# Calling the calculate_circle_area function
radius = 4.5
circle_area = calculate_circle_area(radius)
print("Circle area:", circle_area)
# In Python, a list is a versatile data structure
# that can hold a collection of items. Lists are ordered,
# mutable (meaning you can change their contents), 
# and can contain items of different data types.
# Lists are defined using square brackets []
# and items are separated by commas.
# Here's how you can work with lists:

empty_list = []                     # Create an empty list
numbers = [1, 2, 3, 4, 5]           # Create a list of numbers
fruits = ['apple', 'banana', 'orange']  # Create a list of strings
mixed_list = [1, 'apple', True]      # Create a list with mixed data types

fruits = ['apple', 'banana', 'orange']
print(fruits[0])   # Output: 'apple'
print(fruits[1])   # Output: 'banana'
print(fruits[-1])  # Output: 'orange' (negative indexing starts from the end)

numbers = [1, 2, 3, 4, 5]
subset = numbers[1:4]   # Subset: [2, 3, 4]

fruits = ['apple', 'banana', 'orange']
fruits[1] = 'grape'    # Change 'banana' to 'grape'
fruits.append('kiwi')  # Add 'kiwi' to the end
fruits.insert(1, 'pear')  # Insert 'pear' at index 1
fruits.remove('apple')    # Remove 'apple'

numbers = [3, 1, 4, 1, 5, 9, 2]
numbers.sort()        # Sort the list in ascending order
numbers.reverse()     # Reverse the list
count_1 = numbers.count(1)  # Count occurrences of 1
numbers.pop()         # Remove and return the last element
# List Comprehensions:

# List comprehensions are a concise way to create lists. 
# They allow you to generate new lists by applying an expression
# to each item in an existing iterable.

squares = [x**2 for x in range(1, 6)]  # [1, 4, 9, 16, 25]
# Iterating Over Lists:

fruits = ['apple', 'banana', 'orange']
for fruit in fruits:
    print(fruit)

from collections import deque

# Creating a stack using deque
stack = deque()

# Pushing items onto the stack
stack.append(5)
stack.append("a1")
stack.append(20)

# Peeking at the top item
print("0", stack[0])
print("1", stack[1])
print("2", stack[2])
print("Top item:", stack[-1])

# Popping items from the stack
print("Popped:", stack.pop())
print("Popped:", stack.pop())

# Checking if the stack is empty
print("Is stack empty?", not stack)

# Getting the size of the stack
print("Stack size:", len(stack))
from collections import deque

# Create a deque
deque_example = deque([1, 2, 3, 4, 5], maxlen=5)
print("1.Deque:", deque_example)

# Append and appendleft
deque_example.append(6)
print("2.append:", deque_example)

deque_example.appendleft(0)
print("3.appendleft:", deque_example)

# Extend and extendleft
deque_example.extend([7, 8])
print("4.extend:", deque_example)
deque_example.extendleft([-1, -2])
print("5.extendleft:", deque_example)

# Pop and popleft
popped_item = deque_example.pop()
print("6.pop:", deque_example)

popped_left_item = deque_example.popleft()
print("7.popleft:", deque_example)

# Rotate
deque_example.rotate(1)
print("8.rotate:", deque_example)

# Explanation: Dictionaries in Python are collections 
# of key-value pairs. In this example, a dictionary named user
# is created to store user-related information.
# Values are accessed using the keys.


# Creating a dictionary
user = {
    "username": "johndoe",
    "email": "john@example.com",
    "age": 30
}

# Accessing dictionary values
print("Username:", user["username"])
print("Email:", user["email"])
print("Age:", user["age"])
# Creating a dictionary
student = {
    "name": "Alice",
    "age": 25,
    "courses": ["Math", "Physics"]
}

# Adding a new key-value pair
student["major"] = "Computer Science"

# Removing a key-value pair
del student["age"]

# Checking if a key exists
if "courses" in student:
    print("Student is taking courses:", student["courses"])
# Creating a dictionary
grades = {
    "Math": 90,
    "Science": 85,
    "History": 78
}

# Iterating through dictionary keys and values
for x, y in grades.items():
    print(f"Subject: {x}, Score: {y}")


# Creating a dictionary using dictionary comprehension
numbers_squared = {num: num ** 2 for num in range(1, 6)}

# Printing the squared values
for num, squared in numbers_squared.items():
    print(f"{num} squared is {squared}")
# Explanation: Dictionaries can be nested, allowing
# for more complex data structures. In this example,
#  the employee dictionary contains a nested department
# dictionary.

# These examples demonstrate the usage of dictionaries
# in Python, with explanations tailored for Java experts
# who have experience in web development. 
# They cover creating and accessing dictionaries,
# using dictionary methods, iterating through dictionaries,
# dictionary comprehension, and working with nested dictionaries.






# Creating a nested dictionary
employee = {
    "name": "Alice",
    "department": {
        "name": "Engineering",
        "location": "Building A"
    }
}

# Accessing nested dictionary values
print("Employee:", employee["name"])
print("Department:", employee["department"]["name"])
print("Location:", employee["department"]["location"])
# 1.dict.get(key, default)
# Returns the value for the given key if it exists in
# the dictionary, otherwise returns the specified default value.

grades = {"Math": 90, "Science": 85}
math_grade = grades.get("Math", "N/A")
english_grade = grades.get("English", "N/A")
print("Math Grade:", math_grade)
print("English Grade:", english_grade)

#2. dict.keys()
# Returns a view of the dictionary's keys.

grades = {"Math": 90, "Science": 85}
subject_names = grades.keys()
print("Subjects:", list(subject_names))


#3. dict.values()
# Returns a view of the dictionary's values.

grades = {"Math": 90, "Science": 85}
score_values = grades.values()
print("Scores:", list(score_values))

#4. dict.items()
# Returns a view of the dictionary's key-value pairs as tuples.
grades = {"Math": 90, "Science": 85}
subject_scores = grades.items()
print("Subject Scores:", list(subject_scores))

#5. dict.pop(key, default)
# Removes and returns the value associated with the given key.
#  If the key is not found, it returns the specified default value.

grades = {"Math": 90, "Science": 85}
math_grade = grades.pop("Math")
english_grade = grades.pop("English", -1)
print("Removed Math Grade:", math_grade)
print("Removed English Grade (default):", english_grade)

#6. dict.popitem()
# Removes and returns a random key-value pair as 
# a tuple from the dictionary.

grades = {"Math": 90, "Science": 85, "History": 78}
removed_item = grades.popitem()
print("Removed Item:", removed_item)
print("Updated Dictionary:", grades)

#7. dict.clear()
# dict.clear()
# Removes all items from the dictionary, making it empty.
grades = {"Math": 90, "Science": 85}
grades.clear()
print("Cleared Dictionary:", grades)

# Define a function to calculate the square of a number
def square(x):
    return x ** 2

# List of numbers
numbers = [1, 2, 3, 4, 5]

# Using map to calculate squares of numbers
squared_numbers = map(square, numbers)

# Iterate over the iterator using a loop
for squared_value in squared_numbers:
    print(squared_value)
# Open the file in write mode ('w')
file = open('example.txt', 'w')

# Writing data to the file
file.write('Hello, world!\n')
file.write('This is a sample text.')

# Close the file
file.close()

# The file is now closed and no longer accessible
# Writing to a file with default permissions ('w' mode)
with open('example.txt', 'w') as file:
    file.write('Hello, world!\n')
    file.write('This is a sample text.')

# Reading from a file with default permissions ('r' mode)
with open('example.txt', 'r') as file:
    content = file.read()
    print(content)


# Open the file in write mode ('w')
file = open('example.txt', 'w')

# Writing data to the file
file.write('Hello, world!\n')
file.write('This is a sample text.')

# Close the file
file.close()

# The file is now closed and no longer accessible
import csv
import json

# Writing to a CSV file
data = [['Name', 'Age'], ['Alice', 25], ['Bob', 30]]
with open('data.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(data)

# Reading from a CSV file
with open('data.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        print(row)

# Writing to a JSON file
data = {'name': 'Alice', 'age': 25}
with open('data.json', 'w') as file:
    json.dump(data, file)

# Reading from a JSON file
with open('data.json', 'r') as file:
    data = json.load(file)
    print(data)
try:
    num = int(input("Enter a number: "))
    result = 10 / num
    print(f"Result: {result}")
except ZeroDivisionError:
    print("Cannot divide by zero.")
except ValueError:
    print("Invalid input. Please enter a number.")
except Exception as e:
    print(f"An error occurred: {e}")


# try:
#     print("try")
# except:
#     print("try")# ========raise
# In Python, the raise statement is used to manually raise
# an exception. It allows you to generate and trigger
# exceptions in your code under specific conditions 
# that you define. This can be particularly useful when
# you want to handle specific error cases or enforce
# certain conditions.
# The basic syntax of the raise statement is as follows:


# Define a custom exception class
class MyCustomError(Exception):
    def __init__(self, message):
        self.message = message

try:
    # Attempt to get the user's age from input
    age = int(input("Enter your age: "))
    
    # Check if the age is negative and raise the custom exception if so
    if age < 0:
        raise MyCustomError("Age cannot be negative.")
    
    # If age is not negative, print the age
    print(f"Your age is: {age}")
    
except MyCustomError as e:
    # Handle the custom exception by printing the error message
    print(f"Error: {e}")
    
except ValueError:
    # Handle the case where the input cannot be converted to an integer
    print("Invalid input. Please enter a valid age.")
#new_list = [expression for item in iterable if condition]

squares = [x**2 for x in range(10)]
print(squares)  # Output: [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]

#generator = (expression for item in iterable if condition)

even_generator = (x for x in range(10) if x % 2 == 0)
print(list(even_generator))  # Output: [0, 2, 4, 6, 8]
import numpy as np

# Creating arrays
a = np.array([1, 2, 3])
b = np.array([4, 5, 6])

# Basic operations
c = a + b  # Element-wise addition
d = a * b  # Element-wise multiplication

print("Array a:", a)
print("Array b:", b)
print("Array c (a + b):", c)
print("Array d (a * b):", d)

# Array manipulation
matrix = np.array([[1, 2, 3], [4, 5, 6]])
print("Matrix:")
print(matrix)

# Transposing the matrix
transposed = matrix.T
print("Transposed matrix:")
print(transposed)

# Array statistics
mean_value = np.mean(a)
max_value = np.max(b)
min_index = np.argmin(c)

print("Mean of array a:", mean_value)
print("Maximum value of array b:", max_value)
print("Index of minimum value in array c:", min_index)
import pandas as pd

# Create a DataFrame from a dictionary
data = {'Name': ['Alice', 'Bob', 'Charlie'],
        'Age': [25, 30, 22],
        'City': ['New York', 'Los Angeles', 'Chicago']}
df = pd.DataFrame(data)

# Display the DataFrame
print("DataFrame:")
print(df)

# Accessing columns
names = df['Name']
ages = df['Age']

print("\nNames:")
print(names)
print("\nAges:")
print(ages)

# Filtering data
young_people = df[df['Age'] < 30]
print("\nPeople under 30:")
print(young_people)

# Grouping and aggregation
grouped_cities = df.groupby('City')['Age'].mean()
print("\nAverage age by city:")
print(grouped_cities)
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Create a DataFrame
data = {'Name': ['Alice', 'Bob', 'Charlie'],
        'Age': [25, 30, 22],
        'City': ['New York', 'Los Angeles', 'Chicago']}
df = pd.DataFrame(data)

# Set style for Seaborn
sns.set(style="whitegrid")

# Create a bar plot using Matplotlib
plt.figure(figsize=(8, 6))
plt.bar(df['Name'], df['Age'], color='blue')
plt.xlabel('Name')
plt.ylabel('Age')
plt.title('Age Distribution')
plt.show()

# Create a scatter plot using Seaborn
plt.figure(figsize=(8, 6))
sns.scatterplot(x='Age', y='City', data=df, hue='Name', size='Age', sizes=(50, 200), legend='brief')
plt.xlabel('Age')
plt.ylabel('City')
plt.title('Scatter Plot')
plt.show()
