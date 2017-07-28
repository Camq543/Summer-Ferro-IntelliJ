name := "Ferro Research 2017"

version := "1.0"

scalaVersion := "2.12.2"

classpathTypes += "maven-plugin"

libraryDependencies += "org.nd4j" % "nd4j-native" % "0.8.0" classifier "" classifier "windows-x86_64"

libraryDependencies += "org.nd4j" % "nd4s_2.11" % "0.8.0"
