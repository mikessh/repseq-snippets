def inputFileName = args[0], outputFileName = args[1]

assert args.length == 2

def firstLine = true

new File(outputFileName).withPrintWriter { pwFastq -> 
new File(inputFileName).splitEachLine("\t") { splitLine ->
   if (firstLine) {
       firstLine = false
   } else {
       def seq = splitLine[0], count = splitLine[2]
       if (seq && count)
           pwFastq.println("@" + count  + "\n" + seq + "\n+\n" + ("I" * seq.length()))
   }
}
}