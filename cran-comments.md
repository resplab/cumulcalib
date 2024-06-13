Many thanks for evaluating the package and apologies for the inconveniencs the errors have caused.

Comment: "The Description field is intended to be a (one paragraph) description of what the package does and why it may be useful. Please add more details about the package functionality and implemented methods in your Description text. If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form authors (year) <doi:...> authors (year, ISBN:...) or if those are not available: <https:...> with no space after 'doi:', 'https:' and angle brackets for auto-linking. (If you want to add a title as well please put it in
quotes: "Title")


Response: Description field is now updated, with the publication added. 


Comment: Please write TRUE and FALSE instead of T and F. Please don't use "T" or "F" as vector names.
Response: Done


Comment: Please make sure that you do not change the user's options, par or working directory. If you really have to do so within functions, please ensure with an *immediate* call of on.exit() that the settings are reset when the function is exited.

Response: Done. The oldpar <- par(no.readonly = TRUE) and on.exit(par(oldpar)) are added to the beginning of plot.cumulcalib()


Thanks again for lending your expertise to reviewing our package.

Regards,
Mohsen Sadatsafavi 


## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

