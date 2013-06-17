# Script to Update from Trunk
echo "Please enter your username:"
read username
echo "You are: $username" 

svn up
status=$(svn st |wc -l)

if [ $status -eq 0 ]
then
    svn merge -r408:HEAD svn+ssh://$username@utka42.ethz.ch/rd42/trunk/diamond_analysis
    svn commit -m "Update Branch from Trunk"
    svn up
else
    echo "Please commit before updating"
fi
svn st
