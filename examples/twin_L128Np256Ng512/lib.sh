# Helpers of the step scripts (sourced by config.sh).

# ids of the jobs called NAME that are queued or running with working directory DIR
jobs_in() { squeue -h -u "$USER" -n "$1" -o '%i %Z' 2>/dev/null | awk -v d="$2" '$2 == d {print $1}'; }

# submit_wait NAME DIR SBATCH-ARGS...: submit a job called NAME that runs in DIR, unless one is already queued or
# running there (a rerun after Ctrl-C), then wait for it, checking its state every 30 s. Returns 0 only if the job
# completed.
submit_wait() {
    local name=$1 dir=$2 id state
    shift 2
    id=$(jobs_in "$name" "$dir" | head -1)
    if [ -n "$id" ]; then
        echo "$name: job $id is already queued or running: waiting for it"
    else
        id=$(sbatch --parsable -J "$name" -D "$dir" "$@" | cut -d';' -f1)
        echo "$name: submitted job $id"
    fi
    while :; do
        state=$(sacct -n -X -P -j "$id" -o State 2>/dev/null | head -1)
        case $state in ""|PENDING|RUNNING|REQUEUED|CONFIGURING|COMPLETING|RESIZING|SUSPENDED) sleep 30 ;; *) break ;; esac
    done
    [ "$state" = COMPLETED ] && return 0
    echo "$name: job $id ended $state; see its logs in $dir"
    return 1
}

# trap target for Ctrl-C while waiting: the submitted jobs keep running, and rerunning the step waits for them again
stop_waiting() {
    echo
    echo "stopped waiting; the submitted jobs keep running (squeue -u $USER): rerun this step to wait for them"
    kill $(jobs -p) 2> /dev/null
    exit 130
}
