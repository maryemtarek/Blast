

with open("Query.txt", 'r') as file:
    Query = file.read()

with open('DataBase.txt') as my_file:
    Genes = my_file.readlines()
Genes = [x.strip() for x in Genes]