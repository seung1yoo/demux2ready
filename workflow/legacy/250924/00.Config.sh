PR="$(basename "$PWD")"
path="$(find /bc_storage/BioSeqData/Outsourcing/ -type d -name "$PR" 2>/dev/null || true)"

if [[ -z "$path" ]]; then
    echo -e "\n[ERROR] Project directory not found under /bc_storage/BioSeqData/Outsourcing/: $PR\n" >&2
    exit 1
fi
map=$(find $path -type f -name '*map')

[ -n "$map" ] || { echo -e "\n[ERROR] no files in $map" >&2; exit 1; }
cat $map | tr $'\t' ' ' > map.temp

awk '{$1=$1; print}' map.temp > map.tmp && mv map.tmp map.temp
awk '!seen[$0]++' map.temp > map.tmp && mv map.tmp map.temp

awk 'NR>1 {
    if ($1=="" || $2=="") {
        printf "[ERROR] Empty field at line %d: %s %s\n", NR, $1, $2 > "/dev/stderr"
        err=1
    }
} END {
    if (err) exit 1
}' map.temp || exit 1

# writing config.sh

CONF="config.sh"
echo "PR=\"$PR\"" > $CONF
echo "path=\"$path\"" >> $CONF

outs=$(awk '{print $1}' map.temp | xargs)
reals=$(awk '{print $2}' map.temp | xargs)
rm map.temp
echo "Oname=\"$outs\"" >> $CONF
echo "Rname=\"$reals\"" >> $CONF
echo "SeqType=\"\" # paired-end / single-end" >> $CONF
echo "Mode=\"\" # ngs-only / mb-consortium / dtc-microbiome" >> $CONF
echo "Name=\"\" # Enter customer's name" >> $CONF
echo "Affiliation=\"\" # Enter customer's affiliation" >> $CONF
echo "Date=\"$(date +%y%m%d)\" #YYMMDD ex)250922" >> $CONF
