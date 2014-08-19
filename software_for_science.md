# Let's save ourselves some trouble

For most of us, using good coding practices seems hard to justify. It takes time to learn, time to implement, and hey, our code works just fine. I mean, we're scientists, right?

But the __purpose of software engineering is to save us pain__. It's engineering because we are designing around our own worst tendencies to produce something -- code -- that's robust to the mistakes we are most likely to make.

And yes, doing things this way is going to involve some up-front cost. To save you a world of heartache down the road. I've made all the mistakes, and today I use these techinques in my own work not because I'm a _better_ coder, but because I'm a _humbler_ one: without these habits, I've learned that I'll inevitably produce fragile, unreliable code. 

So let's get started.

## It's about data

For scientists, most of our world revolves around data. How so? Let's look at the activities we tend to spend most of our time on:

- On the basis of old data, we form a hypothesis
- We plan an experiment
    + what data to collect
    + how much to collect
    + power analysis (if you ever do this, seriously, I will buy you lunch)
- We collect data
    + interface with hardware
    + write a program to present stimuli, interact with subjects
- We store data
    + flat files
    + databases
    + video, audio
- We clean data
    + said to be up to 80% of the real work for data scientists
    + get rid of artifacts
    + merging data from various instruments and channels
    + aligning and denoising
    + turning storage formats into variables that are easy to work with
        * matrices
        * data frames
- We analyze data
    + summary statistics
    + hypothesis tests
    + ANOVA/regression
    + quick plots
- We present data
    + fine-tuning figures, statistical tests
- We explain the data
    + write the paper
    + give the talk
    + reproducible research tools

## Scientific programming is special

But look, we're scientists, not programmers. We don't quite have the same sets of goals. In particular

- We don't know where we're going. 
    + This is research. We have to explore a lot of solutions, so no use designing a huge analysis system we'll never use.
    + We do a lot of sloppy, ad hoc analysis
- We are mostly focused on data, not interacting with users.
- We mostly use pre-existing libraries, only occasionally rolling our own
- We generate a lot of plots. Like, _a lot_.
- Our "product" isn't code, it's knowledge.

## Scientific programming is not special

- We need our code to be __correct__ and __robust__
    + unreliable code means unreliable science
- We are just as likely (more so?) than anyone else to write sloppy, hard-to-read, buggy code
- Inevitably, someone else will want to use our code
    + it should be easy to reproduce what we've published
    + it should be easy to understand what it does
    + it should be easy to modify
- __We have to code defensively__

## Why you'll be glad you did

- You will be able to add or remove a single data file/point and instantly re-run all your analyses.
- You will know **exactly** how you normalized the data to produce Figure 2.
- You will be able to **easily** change Figure 2 when Reviewer 3 complains.
- You will be able to return to your code two years later and remember what you did.
- Future rotating graduate students in the lab will not curse your name because they have to modify your code.
- Coding for future projects will become easier over time, since you will be able to confidently reuse old code.
- Good organization will allow you to think at a higher level, so that you spend less time wrestling with your data and more time analyzing.
