while getopts ":ht" opt; do
  case ${opt} in
    h ) # process option h
      ;;
    t ) # process option t
      ;;
    \? ) echo "Usage: cmd [-h] [-t]"
      ;;
  esac
done

echo "$h $t"