# Learning to program

Learning a computer language takes time. Learning how to learn a computer language generally takes much more. Just as with natural languages, the way it typically happens is that you struggle to learn the first new language, struggle a little less with the next, and then somewhere around the third new language, you begin to notice meta-patterns that save you lots of time. From then on, you have a schema for learning, and the features of an individual language become choices among available options, not one-off special cases.

So it's worth our while to think a little about the _generic process_ of writing computer code before we get bogged down in the details of any specific language. The hope is that this information will give you a way to think about learning _any_ computer language, not just the one you end up choosing.

## How to write a computer program

Let's think about what we need to understand to write a computer program, even a small one. For instance, we might want to know:

- What are the inputs to the program?
    + How do I load them (from a file or the internet)?
    + How do I get them from a user?
- What are the outputs of the program?
    + New or modified files?
    + Plots?
    + Statistics?
- What are the steps the program should take to solve the problem (its **algorithms**)?
- What is the information on which the program operates (**data structures**)?
- What problem-solving strategies will we use to combine algorithms and data structures (**design patterns**)?
- Can I break this program into smaller sub-programs, each of which has a single, small job?
- Are there established libraries of functions I can use?
    + What inputs do those libraries require?
    + What data structures do they offer/use?

> The most important skill for programming is the ability to break a large problem into small sub-problems. When problems are small and manageable, it is easier to decide what strategies work best to solve them.

---
## Exercise:

One of the most common ways to analyze spiking neuron data is to create a peri-stimulus or peri-event time histogram (PSTH/PETH). That is, for each event in a set, histogram all spikes in some period of time around that event.

**Describe a program that will construct a PSTH.** Tell what inputs you need and what output you will produce. What steps are necessary? Once you list the steps, try to think through how each might be accomplished. You do not need to write working code. Just a description that another student at your level could understand.

---

## But how do I learn a programming language?

### Read a book
+ Pros: structured introduction, unified perspective
+ Cons: takes a long time, may not cover what you need, not the same as practice

### Study code from the lab
+ Pros: specific to what you want to do, someone you can ask
+ Cons: hard to get started, may teach bad habits

### Some things it helps to know:
+ terminology and philosophy:
    * what's the language "about"
    * hard to Google if you don't know what it's called
+ data types:
    * what kinds of data does the language natively support?
    * what are the "go-to" data types (lists, matrices, data frames, etc.)
+ flow control:
    * how do you branch (`if`/`then`/`else`)?
    * how do you repeat until a test is passed (`do`/`while`)?
    * how do you iterate over a collection (`for`/`in`)?
    * how do you match a with a list of options (`switch`, `case`)?
+ functions:
    * how do I define a reusable block of code?
    * how is information passed to and returned from the block?
    * is a function just another variable type (i.e., a "first-class" function)?
    * can I group functions into modules?
+ objects:
    * what is the object model?
    * what sorts of relationships can exist between objects?
    * what controls which operations can be performed on which objects?
+ i/o:
    * how do I open and close a file?
    * how do I write to a file?
    * how do I get information from users and print to the screen?
+ sugar:
    * special syntax used to perform common operations in programming languages is called "syntactic sugar" (not 100% technically true, but close enough)
    * think "to be" in natural language -- ubiquitous, but always different
    * often implemented with punctuation marks, can be cryptic to new users, hard to Google
    * need to learn enough of this to read basic code

### Google your way to greatness