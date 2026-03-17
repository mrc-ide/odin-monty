knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

knitr::knit_hooks$set(small_margins = function(before, options, envir) {
  if (before) {
    par(mar = c(4, 4, .1, .1))
  }
})

set.seed(42)

lang_output <- function(x, lang) {
  cat(c(sprintf("```%s", lang), x, "```"), sep = "\n")
}

cpp_output <- function(x) {
  lang_output(x, "cpp")
}

r_output <- function(x) {
  lang_output(x, "r")
}

plain_output <- function(x) {
  lang_output(x, "md") # not great, but at least renders nicely
}

model_compile_code <- function(mod_nm, mod_code) {
  c(paste0(mod_nm, " <- odin({"),
    paste0("  ", readLines(mod_code)),
    "})")
}
