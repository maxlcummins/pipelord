while getopts ":io" opt; do
  case ${opt} in
    i ) # process option i
      ;;
    o ) # process option o
      ;;
    \? ) echo "Usage: cmd [-i] [-o]"
      ;;
  esac
done

cat $i

