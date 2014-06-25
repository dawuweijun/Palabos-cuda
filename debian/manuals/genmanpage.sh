#! /bin/sh
#set -x
for name in pvblot mTestDriver vtkEncodeString vtkWrapHierarchy vtkWrapPython vtkWrapPythonInit smTestDriver; do

manname=`echo -n ${name}.1`
cat > $manname <<EOF
.TH NAME 1
.\" NAME $name, SECTION 1
.SH NAME
$name \- $name
.SH SYNOPSIS
.B $name
.br
.SH DESCRIPTION
This manual page documents briefly the
.BR $name
command.

You can find information in the paraview manual page and on the
paraview web page http://www.paraview.org.

.SH AUTHOR
The Debian Scientific Computing Team
EOF
done

