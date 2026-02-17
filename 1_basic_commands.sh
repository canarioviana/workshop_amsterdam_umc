# Connect to the server (Windows PowerShell, Linux terminal, MacOS terminal)
ssh username@ipaddress
# Enter "yes" to save the server identification key to you computer

# Create working directory
mkdir /mnt/4tb_1/workshop_umc/<group_name>/linux_terminal
# Go to working directory
cd /mnt/4tb_1/workshop_umc/<group_name>/linux_terminal

# Print a message in the terminal
echo "Hello World!"

# Create and empty file
touch empty1.txt
touch empty2.txt

# List files
ls

# List files, one per line
ls -1

# List files with properities
ll

# List files with properities, with human readable file sizes
ll -h

# List files with properities and sorted by time
ll -t

# Rename a file
mv empty1.txt Empty1.txt

# Redirect a text to a file
echo "Hello World!" > message.txt # First line
echo "=)"  >> message.txt # Second line

# Read the file content
cat message.txt

# Delete a file
rm message.txt

# Create a directory
mkdir directory

# Create file inside the directory
touch directory/empty3.txt directory/empty4.tsv

# List the directory content
ls directory # all files
ls directory/*.txt # only txt files
ls directory/*.tsv # only tsv files

# Go to directory
cd directory

# Leave the directory (go bach one level)
cd ..

# Move a file to the directory
mv Empty1.txt directory

# List the directory content
ls directory

# Move the file back to current directory
mv directory/Empty1.txt .

# List itens
ls

# Attribute a content to a variable
mymessage="Hello World!"
myfile=directory/empty3.txt

# View the variable content
echo $mymessage
echo $myfile

# Create a loop to present the directory files
for file in directory/*; do
    echo File: $file
done

# Rename a directory
mv directory directory1

# List itens
ls

# Delete the directory
rm -r directory1 # The flag "-r" is required to delete directories

# Go to the main directory
cd /mnt/4tb_1/workshop_umc/