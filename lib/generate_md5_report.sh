#!/bin/bash

DIRECTORY="./"

OUTPUT_FILE="md5_report.txt"

> "$OUTPUT_FILE"

for file in "$DIRECTORY"/*; do
    if [ -f "$file" ]; then
        md5sum "$file" >> "$OUTPUT_FILE"
    fi
done

echo "MD5 해시 리포트가 $OUTPUT_FILE 파일에 생성되었습니다."
