# How to use Processing in IntelliJ Idea? - Useful links
#### [this Stackoverflow question](https://stackoverflow.com/questions/36765288/how-to-use-processing-3-on-intellij-idea)
#### [Written tutorial: Processing in Java](https://happycoding.io/tutorials/java/processing-in-java)

# Steps
1. Download `processing 3.5.4` from [this link](https://processing.org/download).
2. Extract the downloaded archive and add the following libraries to the CLASSPATH:
    - processing-3.5.4/core/library/core.jar
    - processing-3.5.4/core/library/jogl-all.jar
    - processing-3.5.4/core/library/gluegen-rt.jar
    
    NOTE: the following JARs are platform specific, for `Linux x86_64` use:
    - processing-3.5.4/core/library/jogl-all-natives-linux-amd64.jar
    - processing-3.5.4/core/library/gluegen-rt-natives-linux-amd64.jar

    For `Mac OS X` use:
    - processing-3.5.4/core/library/jogl-all-natives-macosx-universal.jar
    - processing-3.5.4/core/library/gluegen-rt-natives-macosx-universal.jar

3. `mkdir statistics frames`
4. Compile
5. Run Main
