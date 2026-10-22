"""Checked subprocess pipes shared by the alignment stages."""

import subprocess
import tempfile
import time


def run_checked_pipe(producer_cmd, consumer_cmd, timeout=None):
    """Run a pipe, checking both exits and terminating children on failure.

    Producer stderr is spooled to disk so a full stderr pipe cannot block
    the producer while the consumer waits for its stdout.
    """
    producer = consumer = None
    deadline = time.monotonic() + timeout if timeout is not None else None

    def remaining():
        return max(0.001, deadline - time.monotonic()) if deadline is not None else None

    with tempfile.TemporaryFile() as producer_stderr:
        try:
            producer = subprocess.Popen(
                producer_cmd, stdout=subprocess.PIPE, stderr=producer_stderr,
            )
            consumer = subprocess.Popen(
                consumer_cmd, stdin=producer.stdout,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            )
            producer.stdout.close()
            _, consumer_stderr = consumer.communicate(timeout=remaining())
            producer.wait(timeout=remaining())
            if producer.returncode:
                producer_stderr.seek(0)
                raise subprocess.CalledProcessError(
                    producer.returncode, producer_cmd,
                    stderr=producer_stderr.read().decode(errors="replace"),
                )
            if consumer.returncode:
                raise subprocess.CalledProcessError(
                    consumer.returncode, consumer_cmd,
                    stderr=consumer_stderr.decode(errors="replace"),
                )
        finally:
            for proc in (consumer, producer):
                if proc is not None:
                    if proc.poll() is None:
                        proc.kill()
                    proc.wait()
            if producer is not None and producer.stdout is not None:
                producer.stdout.close()
