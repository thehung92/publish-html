# publish stuff on the internet

## Let's start simple

Build  a simple plotly html and share it with tiiny.host ~ not working

## host with github

```sh
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote add origin https://github.com/thehung92/publish-html.git
git push -u origin main

# cp the chosen html file to the root directory
cp man/figs/mds3d-merge-clean-cross-impute.html index.html
git add index.html
git commit -m "standalone html"
git push
```

## build sites with rmarkdown

```r
# render the entire site
rmarkdown::render_site()

# render a single file only
rmarkdown::render_site("about.Rmd")
```