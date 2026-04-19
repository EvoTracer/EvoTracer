#!/bin/bash
for file in *.go; do
  if [ "$file" = "transform.sh" ]; then continue; fi
  name="${file%.go}"
  mkdir -p "$name"
  # Change package main to package $name
  sed 's/package main/package '"$name"'/g' "$file" > "$name/$file"
  # Change func main () { to func Main() {
  sed -i 's/func main[[:space:]]*()/func Main()/g' "$name/$file"
  # Remove original file so it doesn't conflict
  rm "$file"
done
