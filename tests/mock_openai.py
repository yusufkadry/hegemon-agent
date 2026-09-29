"""Local stand-in for the OpenAI Responses endpoint (tests only). Replies come from a JSON plan file."""
import json, sys
from http.server import BaseHTTPRequestHandler, HTTPServer

plan_path, log_path, port = sys.argv[1], sys.argv[2], int(sys.argv[3])

class Handler(BaseHTTPRequestHandler):
    def do_POST(self):
        body = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
        with open(log_path, "a") as log:
            log.write(json.dumps({"auth": self.headers.get("Authorization", "")[:12], "body": body}) + "\n")
        plan = json.load(open(plan_path))
        step = plan.pop(0) if len(plan) > 1 else plan[0]
        json.dump(plan, open(plan_path, "w"))
        if step.get("status", 200) != 200:
            payload = {"error": step["error"]}
        else:
            payload = {"status": "completed", "usage": {"input_tokens": 1000, "output_tokens": 200},
                       "output": [{"type": "reasoning"}, {"type": "message", "content": [{"type": "output_text", "text": step["text"]}]}]}
        data = json.dumps(payload).encode()
        self.send_response(step.get("status", 200))
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)
    def log_message(self, *a):
        pass

HTTPServer(("127.0.0.1", port), Handler).serve_forever()
