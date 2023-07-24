import json
import os


def load_expenses():
    if os.path.exists("expenses.json"):
        with open("expenses.json", "r") as file:
            return json.load(file)
    return []


def save_expenses(expenses):
    with open("expenses.json", "w") as file:
        json.dump(expenses, file)

def add_expense(expenses, date, description, amount):
    
    expenses.append({
        'date': date,
        'description': description,
        'amount': amount
    })

def view_expenses(expenses):
    
    print("Expense Tracker")
    print("Date\t\tDescription\tAmount")
    print("--------------------------------------")
    for expense in expenses:
        print(f"{expense['date']}\t{expense['description']}\t${expense['amount']}")

def get_user_input():
   
    while True:
            date = input("Which day of the month is this taken on? (DD): ")
            if 1 <= int(date) <= 31:
                break
            else:
                print("Invalid day. Please enter a day between 1 and 31")        
    description = input("Enter description: ")
    amount = float(input("Enter amount: $"))
    return date, description, amount

def get_income():
    
    income = float(input("Enter your monthly income: $"))
    return income

def view_all_expenses(income, expenses):
    
    view_expenses(expenses)
    disposable_income = calculate_disposable_income(income, expenses)
    print(f"Disposable Income: ${disposable_income:.2f}")
    if disposable_income < 0:
        print("You are spending too much. Slow down!")

def calculate_disposable_income(income, expenses):
    
    total_expenses = sum(expense['amount'] for expense in expenses)
    disposable_income = income - total_expenses
    return disposable_income

def main():
    income = get_income()
    expenses = load_expenses()  
    
    while True:
        print("\n1. Add Expense")
        print("2. View Expenses")
        print("3. View Disposable Income")
        print("4. Compare expenses to income")
        print("5. Exit")
        choice = input("Enter your choice (1/2/3/4/5): ")

        if choice == '1':
            date, description, amount = get_user_input()
            add_expense(expenses, date, description, amount)
        elif choice == '2':
            view_expenses(expenses)
        elif choice == '3':
            disposable_income = calculate_disposable_income(income, expenses)
            print(f"Disposable Income: ${disposable_income:.2f}")
            if disposable_income < 0:
               print("You are spending too much. Slow down!")
        elif choice == '4':
            view_all_expenses(income, expenses)
        elif choice == '5':
            print("Exiting Expense Tracker. Goodbye!")
            break
        else:
            print("Invalid choice. Please try again.")
    
    save_expenses(expenses)  

if __name__ == "__main__":
    main()
