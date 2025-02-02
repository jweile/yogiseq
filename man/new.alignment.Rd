% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/yogiseq.R
\name{new.alignment}
\alias{new.alignment}
\title{creates a new Needleman-Wunsch global alignment}
\usage{
new.alignment(s1, s2)
}
\arguments{
\item{s1}{sequence #1}

\item{s2}{sequence #2}
}
\value{
new alignment with methods:
    <code>getMatrix()</code>
    <code>getDistance()</code>
    <code>getMutations()</code>
    <code>getMappings()</code>
    <code>printAlignment()</code>
    <code>getMatrix()</code>
}
\description{
creates a new Needleman-Wunsch global alignment
}
