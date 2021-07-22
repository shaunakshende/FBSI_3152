for d in ./*/ ; do (cd "$d" && ./Vacuum >/dev/null & disown); done
